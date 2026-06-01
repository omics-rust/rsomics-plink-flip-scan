//! Differential compatibility tests against PLINK 1.9 `--flip-scan`.
//!
//! Each test runs our binary and compares its `.flipscan` output field-by-field
//! (CHR/SNP/A1/A2 strings exact, POS/NEG integers exact, F/R_POS/R_NEG numeric
//! tokens exact — PLINK's `dtoa_g` rounding reproduced — and the NEGSNPS set
//! exact). When the `plink` binary is on PATH we diff live; otherwise we diff
//! against checked-in PLINK 1.9 golden output.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-plink-flip-scan"))
}

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

fn plink_available() -> bool {
    Command::new("plink")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn fields(text: &str) -> Vec<Vec<String>> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(str::to_string).collect())
        .collect()
}

fn assert_fields_equal(ours: &str, reference: &str) {
    let a = fields(ours);
    let b = fields(reference);
    assert_eq!(a.len(), b.len(), "row count differs");
    for (i, (x, y)) in a.iter().zip(&b).enumerate() {
        assert_eq!(x, y, "row {i} differs:\n ours: {x:?}\n ref:  {y:?}");
    }
}

fn run_ours(prefix: &Path) -> String {
    let out = Command::new(ours())
        .arg(prefix)
        .output()
        .expect("run rsomics-plink-flip-scan");
    assert!(
        out.status.success(),
        "rsomics-plink-flip-scan failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("utf8")
}

#[test]
fn small_matches_golden() {
    let ours = run_ours(&golden_dir().join("small"));
    let golden =
        std::fs::read_to_string(golden_dir().join("small.flipscan.golden")).expect("read golden");
    assert_fields_equal(&ours, &golden);
}

#[test]
fn withmiss_matches_golden() {
    let ours = run_ours(&golden_dir().join("withmiss"));
    let golden = std::fs::read_to_string(golden_dir().join("withmiss.flipscan.golden"))
        .expect("read withmiss golden");
    assert_fields_equal(&ours, &golden);
}

#[test]
fn header_is_plink_shape() {
    let ours = run_ours(&golden_dir().join("small"));
    let header: Vec<&str> = ours.lines().next().unwrap().split_whitespace().collect();
    assert_eq!(
        header,
        [
            "CHR", "SNP", "BP", "A1", "A2", "F", "POS", "R_POS", "NEG", "R_NEG", "NEGSNPS"
        ]
    );
}

fn live_diff(name: &str) {
    if !plink_available() {
        eprintln!("plink not on PATH; skipping live differential for {name}");
        return;
    }
    let tmp = tempfile::Builder::new()
        .prefix("plink-flipscan-compat-")
        .tempdir_in(std::env::temp_dir())
        .expect("tempdir");
    let out_prefix = tmp.path().join("ref");
    let status = Command::new("plink")
        .args([
            "--bfile",
            golden_dir().join(name).to_str().unwrap(),
            "--flip-scan",
            "--allow-no-sex",
            "--out",
            out_prefix.to_str().unwrap(),
        ])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("run plink");
    assert!(status.success(), "plink --flip-scan failed");

    let reference =
        std::fs::read_to_string(out_prefix.with_extension("flipscan")).expect("read .flipscan");
    assert_fields_equal(&run_ours(&golden_dir().join(name)), &reference);
}

#[test]
fn small_matches_live_plink() {
    live_diff("small");
}

#[test]
fn withmiss_matches_live_plink() {
    live_diff("withmiss");
}
