use std::path::Path;

use criterion::{Criterion, black_box, criterion_group, criterion_main};
use mmcif_analyze::{
    build_adjacency, build_dense_distance_matrix, build_sparse_bond_adjacency,
    compute_bonds_by_distance, nearest_amino_acids, nearest_atoms,
};
use mmcif_core::{AtomPosition, AtomSite, MmcifStructure, parse_file};

fn load_structure() -> MmcifStructure {
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    let cif_path = manifest.join("../../11AS.cif.gz");
    parse_file(&cif_path).expect("load 11AS")
}

fn anchor_position(structure: &MmcifStructure) -> AtomPosition {
    structure
        .atoms()
        .iter()
        .find_map(|atom| atom.position())
        .expect("structure has coordinates")
}

fn bench_sparse_adjacency(c: &mut Criterion) {
    let structure = load_structure();
    c.bench_function("sparse adjacency 11AS", |b| {
        b.iter(|| {
            let adjacency = build_sparse_bond_adjacency(black_box(&structure), true, 0.35);
            black_box(adjacency);
        });
    });
}

fn bench_dense_adjacency(c: &mut Criterion) {
    let structure = load_structure();
    let natoms = structure.atoms().len();
    c.bench_function("dense adjacency 11AS", |b| {
        b.iter(|| {
            let bonds = compute_bonds_by_distance(black_box(&structure), 0.35);
            let adjacency = build_adjacency(natoms, &bonds);
            black_box(adjacency);
        });
    });
}

fn bench_geometric_bonds(c: &mut Criterion) {
    let structure = load_structure();
    c.bench_function("compute bonds by distance 11AS", |b| {
        b.iter(|| {
            let bonds = compute_bonds_by_distance(black_box(&structure), 0.35);
            black_box(bonds.len());
        });
    });
}

fn bench_dense_distance(c: &mut Criterion) {
    let structure = load_structure();
    let indices: Vec<usize> = structure
        .atoms()
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.position().is_some())
        .map(|(idx, _)| idx)
        .take(128)
        .collect();
    c.bench_function("dense distance (128 atoms)", |b| {
        b.iter(|| {
            let matrix = build_dense_distance_matrix(black_box(&structure), &indices);
            black_box(matrix);
        });
    });
}

fn bench_nearest_atoms(c: &mut Criterion) {
    let structure = load_structure();
    let origin = anchor_position(&structure);
    c.bench_function("nearest atoms search", |b| {
        b.iter(|| {
            let hits = nearest_atoms(black_box(&structure), origin, 5.0, 64);
            black_box(hits);
        });
    });
}

fn bench_nearest_residues(c: &mut Criterion) {
    let structure = load_structure();
    let origin = anchor_position(&structure);
    c.bench_function("nearest residues search", |b| {
        b.iter(|| {
            let hits = nearest_amino_acids(black_box(&structure), origin, 6.0, 64);
            black_box(hits);
        });
    });
}

fn bench_residue_key_prepass(c: &mut Criterion) {
    let structure = load_structure();
    c.bench_function("residue key cache build", |b| {
        b.iter(|| {
            let keys = build_residue_key_cache(structure.atoms());
            black_box(keys.len());
        });
    });
}

fn bench_residue_key_stream(c: &mut Criterion) {
    let structure = load_structure();
    c.bench_function("residue key on-the-fly", |b| {
        b.iter(|| {
            let mut count = 0usize;
            for atom in structure.atoms() {
                if !bench_is_amino_acid(atom.label_comp_id.as_deref()) {
                    continue;
                }
                let _ = derive_residue_key(atom);
                count += 1;
            }
            black_box(count);
        });
    });
}

fn build_residue_key_cache<'a>(atoms: &'a [AtomSite]) -> Vec<Option<BenchResidueKey<'a>>> {
    atoms
        .iter()
        .map(|atom| {
            if bench_is_amino_acid(atom.label_comp_id.as_deref()) {
                Some(derive_residue_key(atom))
            } else {
                None
            }
        })
        .collect()
}

type BenchResidueKey<'a> = (Option<&'a str>, Option<&'a str>, Option<i32>);

fn derive_residue_key<'a>(atom: &'a AtomSite) -> BenchResidueKey<'a> {
    (
        atom.label_comp_id.as_deref(),
        atom.label_asym_id
            .as_deref()
            .or(atom.auth_asym_id.as_deref()),
        atom.label_seq_id.or(atom.auth_seq_id),
    )
}

fn bench_is_amino_acid(comp_id: Option<&str>) -> bool {
    const AMINO_ACIDS: [&str; 21] = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
        "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC",
    ];
    let Some(raw) = comp_id else {
        return false;
    };
    let upper = raw.trim().to_ascii_uppercase();
    AMINO_ACIDS.iter().any(|&aa| aa == upper)
}

criterion_group!(
    analyze_benches,
    bench_sparse_adjacency,
    bench_dense_adjacency,
    bench_geometric_bonds,
    bench_dense_distance,
    bench_nearest_atoms,
    bench_nearest_residues,
    bench_residue_key_prepass,
    bench_residue_key_stream,
);
criterion_main!(analyze_benches);
