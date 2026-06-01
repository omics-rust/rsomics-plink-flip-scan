use clap::Parser;
use rsomics_pgen::Pgen;
use rsomics_plink_flip_scan::{
    DEFAULT_THRESHOLD, DEFAULT_WINDOW, DEFAULT_WINDOW_KB, Params, flip_scan, write_flipscan,
};
use std::fs::File;
use std::io::{self, BufWriter};
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "rsomics-plink-flip-scan",
    about = "PLINK1 --flip-scan: LD-sign strand-inconsistency QC scan",
    version
)]
struct Cli {
    /// Path prefix for the .bed/.bim/.fam fileset (without extension).
    bfile: PathBuf,

    /// Max variant-count distance scanned each side (plink --flip-scan-window).
    #[arg(long = "flip-scan-window", default_value_t = DEFAULT_WINDOW)]
    window: usize,

    /// Max kb distance scanned (plink --flip-scan-window-kb).
    #[arg(long = "flip-scan-window-kb", default_value_t = DEFAULT_WINDOW_KB)]
    window_kb: f64,

    /// Min correlation for a partner to count (plink --flip-scan-threshold).
    #[arg(long = "flip-scan-threshold", default_value_t = DEFAULT_THRESHOLD)]
    threshold: f64,

    /// Write the report to <OUT>.flipscan instead of stdout (plink --out).
    #[arg(short = 'o', long)]
    out: Option<PathBuf>,
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: Cli) -> anyhow::Result<()> {
    let pgen = Pgen::load(&cli.bfile)?;
    let params = Params {
        window: cli.window,
        window_kb: cli.window_kb,
        threshold: cli.threshold,
    };
    let records = flip_scan(&pgen, &params);

    match cli.out {
        Some(prefix) => {
            let path = prefix.with_extension("flipscan");
            let mut w = BufWriter::new(File::create(path)?);
            write_flipscan(&records, &mut w)?;
        }
        None => {
            let stdout = io::stdout();
            let mut w = BufWriter::new(stdout.lock());
            write_flipscan(&records, &mut w)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_definition_is_valid() {
        <Cli as clap::CommandFactory>::command().debug_assert();
    }
}
