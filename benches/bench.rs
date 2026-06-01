use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_pgen::Pgen;
use rsomics_plink_flip_scan::{Params, flip_scan};
use std::hint::black_box;
use std::path::PathBuf;

fn bench_flip_scan(c: &mut Criterion) {
    let prefix = std::env::var("FLIPSCAN_BENCH_BFILE")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/small"));
    let pgen = Pgen::load(&prefix).expect("load fileset");
    let params = Params::default();
    c.bench_function("flip_scan", |b| {
        b.iter(|| flip_scan(black_box(&pgen), black_box(&params)));
    });
}

criterion_group!(benches, bench_flip_scan);
criterion_main!(benches);
