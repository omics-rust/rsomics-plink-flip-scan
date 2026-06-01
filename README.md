# rsomics-plink-flip-scan

PLINK1 `--flip-scan`: an LD-sign strand-inconsistency QC scan. For each SNP it
compares the genotype correlation (r) with nearby SNPs between cases and
controls; a partner whose correlation flips sign between the two groups
(negative where the rest of the block is positive) is flagged as a likely
strand-flip or allele-coding error. Output is plink's `.flipscan` table.

```
rsomics-plink-flip-scan PREFIX [--flip-scan-window N] [--flip-scan-window-kb X] \
                        [--flip-scan-threshold X] [-o OUT]
```

`PREFIX` is the `.bed`/`.bim`/`.fam` fileset prefix. Without `-o`, the report
goes to stdout; with `-o OUT` it is written to `OUT.flipscan`.

Columns: `CHR SNP BP A1 A2 F POS R_POS NEG R_NEG NEGSNPS`, where `POS`/`NEG`
count the in-window partners whose case/control correlations agree/disagree in
sign (gated on `max(|r_case|, |r_ctrl|) >= threshold`), `R_POS`/`R_NEG` are the
mean of `(|r_case| + |r_ctrl|) / 2` over those partners, and `NEGSNPS` lists the
negative-match SNPs.

Defaults match plink: window 10 variants each side bounded by 1000 kb,
correlation threshold 0.5.

## Performance

Byte-for-byte plink-identical output, faster on the hot path: on a 50 000-SNP ×
2 000-sample fileset (PLINK v1.9.0-b.7.7, Xeon Gold 6238R, single core), CPU
time `0.61 s` vs plink `0.65 s` (1.07×), via packed-genotype bitplane LD
(four AND-popcounts per 64-sample word) and a streaming window that keeps only
the active SNPs resident.

## Origin

This crate is an independent Rust reimplementation of PLINK1 `--flip-scan`
based on:
- The published method (Chang et al. 2015, PLINK 1.9, doi:10.1186/s13742-015-0047-8;
  Purcell et al. 2007, PLINK 1, doi:10.1086/519795)
- The public PLINK 1.9 LD / file-format specs
  (https://www.cog-genomics.org/plink/1.9/ld , https://www.cog-genomics.org/plink/1.9/formats)
- Black-box behaviour testing against the upstream plink binary

No source code from the GPL upstream was used as reference during
implementation. Test fixtures are independently generated.

License: MIT OR Apache-2.0.
Upstream credit: PLINK 1.9 (Christopher Chang et al., GPLv3).
