use std::path::Path;

use criterion::{Criterion, black_box, criterion_group, criterion_main};
use mmcif_core::parse_file;

fn bench_parse_11as(c: &mut Criterion) {
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    let cif_path = manifest.join("../../11AS.cif.gz");
    c.bench_function("mmcif_core parse 11AS", |b| {
        b.iter(|| {
            let structure = parse_file(&cif_path).expect("parse cif");
            black_box(structure.atom_count());
        });
    });
}

criterion_group!(core_benches, bench_parse_11as);
criterion_main!(core_benches);
