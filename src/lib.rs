//! PLINK1 `--flip-scan`: LD-sign strand-inconsistency QC.
//!
//! For each index SNP it scans nearby SNPs (within a variant-count and bp
//! window on the same chromosome) and, for every partner, computes the
//! genotype correlation r separately in cases (`R_A`) and controls (`R_U`).
//! A partner is a *match* when `max(|R_A|, |R_U|) >= threshold`; a match with
//! agreeing signs is positive, with opposing signs negative. A negative match
//! flags a likely strand flip or allele-coding error between the two groups.
//!
//! The reported `R_POS`/`R_NEG` are the mean of `(|R_A| + |R_U|) / 2` over the
//! matched partners.

#![allow(clippy::cast_precision_loss)]

use rsomics_pgen::Pgen;
use std::fmt::Write as FmtWrite;
use std::io::{self, Write};

pub const DEFAULT_WINDOW: usize = 10;
pub const DEFAULT_WINDOW_KB: f64 = 1000.0;
pub const DEFAULT_THRESHOLD: f64 = 0.5;

pub struct Params {
    pub window: usize,
    pub window_kb: f64,
    pub threshold: f64,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            window: DEFAULT_WINDOW,
            window_kb: DEFAULT_WINDOW_KB,
            threshold: DEFAULT_THRESHOLD,
        }
    }
}

pub struct FlipScanRecord {
    pub chrom: String,
    pub snp: String,
    pub bp: u64,
    pub a1: String,
    pub a2: String,
    pub maf: f64,
    pub pos: u32,
    pub r_pos: f64,
    pub neg: u32,
    pub r_neg: f64,
    pub negsnps: Vec<String>,
}

/// Per-byte (4 samples) decode tables mapping a packed .bed byte to the 4-bit
/// nibble of each plane: `p1` = dosage ≥ 1 (HomA1/Het), `p2` = dosage 2
/// (HomA1), `present` = non-missing.
struct DecodeLut {
    p1: [u8; 256],
    p2: [u8; 256],
    present: [u8; 256],
}

impl DecodeLut {
    fn build() -> Self {
        let mut p1 = [0u8; 256];
        let mut p2 = [0u8; 256];
        let mut present = [0u8; 256];
        for byte in 0usize..256 {
            for lane in 0..4 {
                match (byte >> (lane * 2)) & 0b11 {
                    0b00 => {
                        p1[byte] |= 1 << lane;
                        p2[byte] |= 1 << lane;
                        present[byte] |= 1 << lane;
                    }
                    0b10 => {
                        p1[byte] |= 1 << lane;
                        present[byte] |= 1 << lane;
                    }
                    0b11 => present[byte] |= 1 << lane,
                    _ => {}
                }
            }
        }
        Self { p1, p2, present }
    }
}

/// One SNP's genotypes in a phenotype group, compacted to that group's samples.
///
/// Dosage `d ∈ {0,1,2}` splits into `p1 = (d >= 1)`, `p2 = (d == 2)`; `present`
/// marks non-missing. Correlation in the group is then AND-popcounts over 64
/// samples per word — the plink LD trick — with per-SNP `sum`/`sumsq` cached.
#[derive(Clone, Default)]
struct SnpPlane {
    p1: Vec<u64>,
    p2: Vec<u64>,
    present: Vec<u64>,
    sum: f64,
    sumsq: f64,
    has_missing: bool,
}

impl SnpPlane {
    fn alloc(words: usize) -> Self {
        Self {
            p1: vec![0; words],
            p2: vec![0; words],
            present: vec![0; words],
            sum: 0.0,
            sumsq: 0.0,
            has_missing: false,
        }
    }
}

/// Streaming flip-scan engine: decodes each SNP into a ring buffer of the last
/// `max_dist + 1` SNPs (one set of compacted planes per phenotype group), so
/// only the active window stays resident — never the whole genotype matrix.
struct Scan {
    words_all: usize,
    group_words: [usize; 2],
    group_n: [u32; 2],
    plans: [CompactPlan; 2],
    founder_mask: Vec<u64>,
    ring: usize,
    /// `[group][slot]` compacted planes; slot = snp index mod `ring`.
    slots: [Vec<SnpPlane>; 2],
    full_p1: Vec<u64>,
    full_p2: Vec<u64>,
    full_present: Vec<u64>,
    lut: DecodeLut,
    pub freq: Vec<f64>,
}

impl Scan {
    fn new(
        pgen: &Pgen,
        case_mask: Vec<u64>,
        ctrl_mask: Vec<u64>,
        founder_mask: &[u64],
        max_dist: usize,
    ) -> Self {
        let n_samp = pgen.n_samples();
        let words_all = n_samp.div_ceil(64);
        let group_n = [
            case_mask.iter().map(|w| w.count_ones()).sum(),
            ctrl_mask.iter().map(|w| w.count_ones()).sum(),
        ];
        let plans = [CompactPlan::new(&case_mask), CompactPlan::new(&ctrl_mask)];
        let group_words = [plans[0].out_words, plans[1].out_words];
        let ring = max_dist + 1;
        let slots = [
            (0..ring).map(|_| SnpPlane::alloc(group_words[0])).collect(),
            (0..ring).map(|_| SnpPlane::alloc(group_words[1])).collect(),
        ];
        Self {
            words_all,
            group_words,
            group_n,
            plans,
            founder_mask: founder_mask.to_vec(),
            ring,
            slots,
            full_p1: vec![0; words_all],
            full_p2: vec![0; words_all],
            full_present: vec![0; words_all],
            lut: DecodeLut::build(),
            freq: vec![0.0; pgen.n_variants()],
        }
    }

    /// Decode SNP `i` into the full-sample scratch planes, record its founder
    /// allele frequency, then compact into both groups' ring slots.
    fn load(&mut self, pgen: &Pgen, i: usize) {
        let n_samp = pgen.n_samples();
        let row = pgen.variant_row(i);
        self.full_p1.iter_mut().for_each(|w| *w = 0);
        self.full_p2.iter_mut().for_each(|w| *w = 0);
        self.full_present.iter_mut().for_each(|w| *w = 0);
        for (b, &byte) in row.iter().enumerate() {
            let byte = byte as usize;
            let (w, sh) = (b / 16, (b % 16) * 4);
            self.full_p1[w] |= u64::from(self.lut.p1[byte]) << sh;
            self.full_p2[w] |= u64::from(self.lut.p2[byte]) << sh;
            self.full_present[w] |= u64::from(self.lut.present[byte]) << sh;
        }
        let tail = n_samp % 64;
        if tail != 0 {
            let mask = (1u64 << tail) - 1;
            let last = self.words_all - 1;
            self.full_p1[last] &= mask;
            self.full_p2[last] &= mask;
            self.full_present[last] &= mask;
        }

        let (mut a1, mut nobs, mut all_present) = (0u64, 0u64, 0u64);
        for k in 0..self.words_all {
            let m = self.founder_mask[k];
            a1 += (self.full_p1[k] & m).count_ones() as u64
                + (self.full_p2[k] & m).count_ones() as u64;
            nobs += (self.full_present[k] & m).count_ones() as u64;
            all_present += self.full_present[k].count_ones() as u64;
        }
        self.freq[i] = if nobs == 0 {
            0.0
        } else {
            a1 as f64 / (2.0 * nobs as f64)
        };
        let any_missing = all_present != n_samp as u64;

        let slot = i % self.ring;
        for g in 0..2 {
            let gw = self.group_words[g];
            let pl = &mut self.slots[g][slot];
            compact_with(&self.full_p1, &self.plans[g], &mut pl.p1);
            compact_with(&self.full_p2, &self.plans[g], &mut pl.p2);
            let (mut s, mut ss) = (0u64, 0u64);
            for k in 0..gw {
                let c1 = pl.p1[k].count_ones() as u64;
                let c2 = pl.p2[k].count_ones() as u64;
                s += c1 + c2;
                ss += c1 + 3 * c2;
            }
            pl.sum = s as f64;
            pl.sumsq = ss as f64;
            // present is only consulted on the missing path; compact it lazily.
            pl.has_missing = if any_missing {
                compact_with(&self.full_present, &self.plans[g], &mut pl.present);
                pl.present
                    .iter()
                    .map(|w| w.count_ones() as u64)
                    .sum::<u64>()
                    != u64::from(self.group_n[g])
            } else {
                false
            };
        }
    }

    #[inline]
    fn corr_case(&self, j: usize, i: usize) -> Option<f64> {
        self.corr(0, j, i)
    }
    #[inline]
    fn corr_ctrl(&self, j: usize, i: usize) -> Option<f64> {
        self.corr(1, j, i)
    }

    #[inline]
    fn corr(&self, g: usize, j: usize, i: usize) -> Option<f64> {
        let a = &self.slots[g][j % self.ring];
        let b = &self.slots[g][i % self.ring];
        let w = self.group_words[g];
        if !(a.has_missing || b.has_missing) {
            let mut sxy = 0u64;
            for k in 0..w {
                let (x1, x2) = (a.p1[k], a.p2[k]);
                sxy += (x1 & b.p1[k]).count_ones() as u64
                    + (x1 & b.p2[k]).count_ones() as u64
                    + (x2 & b.p1[k]).count_ones() as u64
                    + (x2 & b.p2[k]).count_ones() as u64;
            }
            return finish_r(
                f64::from(self.group_n[g]),
                a.sum,
                b.sum,
                a.sumsq,
                b.sumsq,
                sxy as f64,
            );
        }
        let (mut n, mut sx, mut sy, mut sxx, mut syy, mut sxy) =
            (0u64, 0u64, 0u64, 0u64, 0u64, 0u64);
        for k in 0..w {
            let both = a.present[k] & b.present[k];
            n += both.count_ones() as u64;
            let (a1, a2) = (a.p1[k] & both, a.p2[k] & both);
            let (b1, b2) = (b.p1[k] & both, b.p2[k] & both);
            sx += a1.count_ones() as u64 + a2.count_ones() as u64;
            sy += b1.count_ones() as u64 + b2.count_ones() as u64;
            sxx += a1.count_ones() as u64 + 3 * a2.count_ones() as u64;
            syy += b1.count_ones() as u64 + 3 * b2.count_ones() as u64;
            sxy += (a1 & b1).count_ones() as u64
                + (a1 & b2).count_ones() as u64
                + (a2 & b1).count_ones() as u64
                + (a2 & b2).count_ones() as u64;
        }
        finish_r(
            n as f64, sx as f64, sy as f64, sxx as f64, syy as f64, sxy as f64,
        )
    }
}

/// Precomputed compaction plan for one group: the source words that carry any
/// selected bit, with each word's mask and bit count, so `compact_with` skips
/// empty words and never recomputes the popcount.
struct CompactPlan {
    /// `(source_word_index, mask_word, bit_count)`.
    words: Vec<(usize, u64, u32)>,
    out_words: usize,
}

impl CompactPlan {
    fn new(mask: &[u64]) -> Self {
        let words: Vec<(usize, u64, u32)> = mask
            .iter()
            .enumerate()
            .filter(|&(_, &m)| m != 0)
            .map(|(i, &m)| (i, m, m.count_ones()))
            .collect();
        let total: u32 = words.iter().map(|&(_, _, c)| c).sum();
        Self {
            words,
            out_words: (total as usize).div_ceil(64).max(1),
        }
    }
}

/// Extract the plan's selected bits of `src` and pack them densely into `dst`.
fn compact_with(src: &[u64], plan: &CompactPlan, dst: &mut [u64]) {
    dst.iter_mut().for_each(|w| *w = 0);
    let mut acc = 0u64;
    let mut filled = 0u32;
    let mut out = 0usize;
    for &(idx, m, cnt) in &plan.words {
        let bits = pext(src[idx], m);
        acc |= bits << filled;
        if filled + cnt >= 64 {
            dst[out] = acc;
            out += 1;
            acc = if filled == 0 {
                0
            } else {
                bits >> (64 - filled)
            };
            filled = (filled + cnt) - 64;
        } else {
            filled += cnt;
        }
    }
    if out < dst.len() {
        dst[out] = acc;
    }
}

/// Parallel bit extract: gather the bits of `src` selected by `mask` into the
/// low `popcount(mask)` bits. Uses the BMI2 `pext` instruction where present,
/// with a portable bit-by-bit fallback elsewhere.
#[inline]
fn pext(src: u64, mask: u64) -> u64 {
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("bmi2") {
            return unsafe { core::arch::x86_64::_pext_u64(src, mask) };
        }
    }
    let mut res = 0u64;
    let mut m = mask;
    let mut out = 0;
    while m != 0 {
        let bit = m & m.wrapping_neg();
        if src & bit != 0 {
            res |= 1 << out;
        }
        out += 1;
        m &= m - 1;
    }
    res
}

#[inline]
fn finish_r(nf: f64, sx: f64, sy: f64, sxx: f64, syy: f64, sxy: f64) -> Option<f64> {
    if nf < 2.0 {
        return None;
    }
    let cov = sxy - sx * sy / nf;
    let vx = sxx - sx * sx / nf;
    let vy = syy - sy * sy / nf;
    if vx <= 0.0 || vy <= 0.0 {
        return None;
    }
    Some(cov / (vx * vy).sqrt())
}

/// Build a `words`-sized sample mask selecting the indices for which `keep` is
/// true.
fn sample_mask(keep: impl Iterator<Item = bool>, words: usize) -> Vec<u64> {
    let mut mask = vec![0u64; words];
    for (s, k) in keep.enumerate() {
        if k {
            mask[s / 64] |= 1u64 << (s % 64);
        }
    }
    mask
}

/// Per-SNP correlation tally, accumulated symmetrically as each unordered pair
/// is scored once. Partners are recorded in ascending SNP order because the
/// scan visits pairs `(j, i)` with `j < i` in increasing `i`.
#[derive(Default)]
struct Tally {
    pos_sum: f64,
    pos_ct: u32,
    neg_sum: f64,
    neg_ct: u32,
    negsnps: Vec<u32>,
}

#[must_use]
pub fn flip_scan(pgen: &Pgen, p: &Params) -> Vec<FlipScanRecord> {
    let n = pgen.n_variants();
    let max_dist = p.window.saturating_sub(1);
    let kb_span = (p.window_kb * 1000.0) as u64;

    let n_samp = pgen.n_samples();
    let words = n_samp.div_ceil(64);
    let case_mask = sample_mask(pgen.samples.iter().map(|s| s.phen == "2"), words);
    let ctrl_mask = sample_mask(pgen.samples.iter().map(|s| s.phen == "1"), words);
    let founder_mask = sample_mask(
        pgen.samples.iter().map(|s| s.pid == "0" && s.mid == "0"),
        words,
    );

    let mut scan = Scan::new(pgen, case_mask, ctrl_mask, &founder_mask, max_dist);
    let mut tally: Vec<Tally> = (0..n).map(|_| Tally::default()).collect();

    for i in 0..n {
        scan.load(pgen, i);
        let lo = i.saturating_sub(max_dist);
        for j in lo..i {
            let vi = &pgen.variants[i];
            let vj = &pgen.variants[j];
            if vj.chrom != vi.chrom || vi.pos.abs_diff(vj.pos) > kb_span {
                continue;
            }
            let (Some(ra), Some(ru)) = (scan.corr_case(j, i), scan.corr_ctrl(j, i)) else {
                continue;
            };
            if ra.abs().max(ru.abs()) < p.threshold {
                continue;
            }
            let mean_abs = (ra.abs() + ru.abs()) / 2.0;
            if ra.signum() == ru.signum() {
                tally[i].pos_sum += mean_abs;
                tally[i].pos_ct += 1;
                tally[j].pos_sum += mean_abs;
                tally[j].pos_ct += 1;
            } else {
                tally[i].neg_sum += mean_abs;
                tally[i].neg_ct += 1;
                tally[i].negsnps.push(j as u32);
                tally[j].neg_sum += mean_abs;
                tally[j].neg_ct += 1;
                tally[j].negsnps.push(i as u32);
            }
        }
    }

    tally
        .into_iter()
        .enumerate()
        .map(|(i, t)| {
            let vi = &pgen.variants[i];
            FlipScanRecord {
                chrom: vi.chrom.clone(),
                snp: vi.id.clone(),
                bp: vi.pos,
                a1: vi.a1.clone(),
                a2: vi.a2.clone(),
                maf: scan.freq[i],
                pos: t.pos_ct,
                r_pos: if t.pos_ct == 0 {
                    f64::NAN
                } else {
                    t.pos_sum / f64::from(t.pos_ct)
                },
                neg: t.neg_ct,
                r_neg: if t.neg_ct == 0 {
                    f64::NAN
                } else {
                    t.neg_sum / f64::from(t.neg_ct)
                },
                negsnps: t
                    .negsnps
                    .iter()
                    .map(|&j| pgen.variants[j as usize].id.clone())
                    .collect(),
            }
        })
        .collect()
}

const CHR_W: usize = 6;
const BP_W: usize = 13;
const A_W: usize = 5;
const F_W: usize = 9;
const CT_W: usize = 7;
const R_W: usize = 9;

/// SNP column field width (right-justified, includes its left separator),
/// matching plink's report convention for the marker-ID column.
fn snp_width(max_id: usize) -> usize {
    if max_id <= 4 { 5 } else { max_id + 3 }
}

/// Write the records in plink's `.flipscan` text layout (right-justified
/// columns, `NA` for empty correlations, `|`-joined negative-match SNP list).
pub fn write_flipscan<W: Write>(records: &[FlipScanRecord], out: &mut W) -> io::Result<()> {
    let sw = snp_width(records.iter().map(|r| r.snp.len()).max().unwrap_or(0));
    let mut line = String::with_capacity(128);
    for (label, width) in [
        ("CHR", CHR_W),
        ("SNP", sw),
        ("BP", BP_W),
        ("A1", A_W),
        ("A2", A_W),
        ("F", F_W),
        ("POS", CT_W),
        ("R_POS", R_W),
        ("NEG", CT_W),
        ("R_NEG", R_W),
    ] {
        pad_str(&mut line, label, width);
    }
    line.push_str(" NEGSNPS\n");
    out.write_all(line.as_bytes())?;
    let mut num = String::with_capacity(16);
    for r in records {
        line.clear();
        pad_str(&mut line, &r.chrom, CHR_W);
        pad_str(&mut line, &r.snp, sw);
        num.clear();
        write!(num, "{}", r.bp).unwrap();
        pad_str(&mut line, &num, BP_W);
        pad_str(&mut line, &r.a1, A_W);
        pad_str(&mut line, &r.a2, A_W);
        num.clear();
        fmt_g_into(&mut num, r.maf);
        pad_str(&mut line, &num, F_W);
        num.clear();
        write!(num, "{}", r.pos).unwrap();
        pad_str(&mut line, &num, CT_W);
        num.clear();
        fmt_na_into(&mut num, r.r_pos);
        pad_str(&mut line, &num, R_W);
        num.clear();
        write!(num, "{}", r.neg).unwrap();
        pad_str(&mut line, &num, CT_W);
        num.clear();
        fmt_na_into(&mut num, r.r_neg);
        pad_str(&mut line, &num, R_W);
        line.push(' ');
        for (k, s) in r.negsnps.iter().enumerate() {
            if k != 0 {
                line.push('|');
            }
            line.push_str(s);
        }
        line.push('\n');
        out.write_all(line.as_bytes())?;
    }
    Ok(())
}

/// Right-justify `s` into `width` (space-padded) and append to `line`.
fn pad_str(line: &mut String, s: &str, width: usize) {
    for _ in s.len()..width {
        line.push(' ');
    }
    line.push_str(s);
}

/// A fixed stack buffer implementing `fmt::Write`, for formatting a shortest
/// decimal without heap allocation.
struct StackBuf {
    b: [u8; 32],
    len: usize,
}

impl FmtWrite for StackBuf {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        let bytes = s.as_bytes();
        self.b[self.len..self.len + bytes.len()].copy_from_slice(bytes);
        self.len += bytes.len();
        Ok(())
    }
}

/// plink's `dtoa_g` at 3 significant figures, written into `buf`.
///
/// plink does not round the exact f64: it first takes the shortest round-trip
/// decimal (Rust's default `Display`), then rounds *that* decimal to 3 sig figs
/// with round-half-to-even. So 0.2875 (whose true f64 is 0.28749…) rounds to
/// 0.288, not the correctly-rounded 0.287 — that's what makes the F / R columns
/// byte-identical. Scientific when the exponent is `< -4` or `>= 3`.
fn fmt_g_into(buf: &mut String, x: f64) {
    if x == 0.0 {
        buf.push('0');
        return;
    }
    if x.is_nan() {
        buf.push_str("nan");
        return;
    }
    let neg = x < 0.0;
    let mut sb = StackBuf { b: [0; 32], len: 0 };
    write!(sb, "{}", x.abs()).unwrap();
    let shortest = std::str::from_utf8(&sb.b[..sb.len]).unwrap();

    const P: i32 = 3;
    let (digits, exp10) = decimal_digits(shortest);
    let (rounded, carry_exp) = round_half_even(&digits, P as usize);
    let exp = exp10 + carry_exp;

    if neg {
        buf.push('-');
    }
    if !(-4..P).contains(&exp) {
        let frac = rounded[1..].trim_end_matches('0');
        buf.push_str(&rounded[..1]);
        if !frac.is_empty() {
            buf.push('.');
            buf.push_str(frac);
        }
        buf.push('e');
        buf.push(if exp < 0 { '-' } else { '+' });
        let e = exp.unsigned_abs();
        if e < 10 {
            buf.push('0');
        }
        write!(buf, "{e}").unwrap();
    } else if exp < 0 {
        buf.push_str("0.");
        for _ in 0..(-exp - 1) {
            buf.push('0');
        }
        let t = rounded.trim_end_matches('0');
        buf.push_str(if t.is_empty() { "0" } else { t });
    } else {
        let intlen = (exp + 1) as usize;
        if rounded.len() <= intlen {
            buf.push_str(&rounded);
            for _ in 0..(intlen - rounded.len()) {
                buf.push('0');
            }
        } else {
            buf.push_str(&rounded[..intlen]);
            let frac = rounded[intlen..].trim_end_matches('0');
            if !frac.is_empty() {
                buf.push('.');
                buf.push_str(frac);
            }
        }
    }
}

fn fmt_na_into(buf: &mut String, x: f64) {
    if x.is_nan() {
        buf.push_str("NA");
    } else {
        fmt_g_into(buf, x);
    }
}

/// Decompose a shortest-decimal string into its significant digit string (no
/// leading zeros, no point) and the base-10 exponent of the first digit.
fn decimal_digits(s: &str) -> (String, i32) {
    let (int_part, frac_part) = s.split_once('.').unwrap_or((s, ""));
    let raw: String = format!("{int_part}{frac_part}");
    let lead = raw.chars().take_while(|&c| c == '0').count();
    let digits: String = raw[lead..].trim_end_matches('0').to_string();
    let digits = if digits.is_empty() {
        "0".into()
    } else {
        digits
    };
    let exp = if int_part.trim_start_matches('0').is_empty() {
        -(1 + frac_part.chars().take_while(|&c| c == '0').count() as i32)
    } else {
        int_part.trim_start_matches('0').len() as i32 - 1
    };
    (digits, exp)
}

/// Round a significant-digit string to `p` digits, half-to-even. Returns the
/// rounded digit string and an exponent shift (1 if a carry grew a digit).
fn round_half_even(digits: &str, p: usize) -> (String, i32) {
    let bytes = digits.as_bytes();
    if bytes.len() <= p {
        return (digits.to_string(), 0);
    }
    let mut kept: Vec<u8> = bytes[..p].iter().map(|b| b - b'0').collect();
    let next = bytes[p] - b'0';
    let rest_nonzero = bytes[p + 1..].iter().any(|&b| b != b'0');
    let round_up = next > 5 || (next == 5 && (rest_nonzero || (kept[p - 1] & 1) == 1));
    let mut carry_exp = 0;
    if round_up {
        let mut i = p;
        loop {
            if i == 0 {
                kept.insert(0, 1);
                kept.pop();
                carry_exp = 1;
                break;
            }
            i -= 1;
            if kept[i] == 9 {
                kept[i] = 0;
            } else {
                kept[i] += 1;
                break;
            }
        }
    }
    let s: String = kept.iter().map(|d| (d + b'0') as char).collect();
    (s, carry_exp)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g(x: f64) -> String {
        let mut s = String::new();
        fmt_g_into(&mut s, x);
        s
    }

    fn na(x: f64) -> String {
        let mut s = String::new();
        fmt_na_into(&mut s, x);
        s
    }

    #[test]
    fn g_matches_plink() {
        assert_eq!(g(0.29), "0.29");
        assert_eq!(g(0.577543), "0.578");
        assert_eq!(g(0.531208), "0.531");
        assert_eq!(g(0.59033), "0.59");
        assert_eq!(g(0.700219), "0.7");
        assert_eq!(g(0.3), "0.3");
        assert_eq!(g(0.0), "0");
        assert_eq!(na(f64::NAN), "NA");
    }

    /// plink rounds the shortest decimal, not the exact f64, half-to-even:
    /// 0.2875 → 0.288, 0.3225 → 0.322, even though the true f64s sit on the
    /// other side of the half.
    #[test]
    fn g_rounds_shortest_decimal_half_even() {
        assert_eq!(g(115.0 / 400.0), "0.288"); // 0.2875
        assert_eq!(g(123.0 / 400.0), "0.308"); // 0.3075
        assert_eq!(g(129.0 / 400.0), "0.322"); // 0.3225
        assert_eq!(g(119.0 / 400.0), "0.298"); // 0.2975
        assert_eq!(g(109.0 / 400.0), "0.272"); // 0.2725
    }

    #[test]
    fn g_scientific_and_carry() {
        assert_eq!(g(0.0001234), "0.000123");
        assert_eq!(g(0.00001234), "1.23e-05");
        assert_eq!(g(999.9), "1e+03");
        assert_eq!(g(12.34), "12.3");
        assert_eq!(g(1.0), "1");
    }
}
