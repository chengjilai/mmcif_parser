#![forbid(unsafe_code)]
//! Spatial search helpers built on top of `mmcif_core` data structures.
//! The API is intentionally small: feed in a `MmcifStructure`, a query point,
//! a distance cutoff, and a hit cap to retrieve the closest atoms or amino-acid
//! residues.

use std::cmp::Ordering;
use std::collections::HashMap;

use mmcif_core::atom::{BondType, Hybridization};
use mmcif_core::{AtomPosition, AtomSite, BondEdge, MmcifStructure, UnitCell};

/// Result entry for [`nearest_atoms`].
#[derive(Debug, Clone)]
pub struct AtomHit<'a> {
    pub index: usize,
    pub distance: f64,
    pub atom: &'a AtomSite,
}

/// Resolved bond mapped to atom indices.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedBond {
    pub first: usize,
    pub second: usize,
    pub bond_type: Option<BondType>,
}

/// Result entry for [`nearest_amino_acids`].
#[derive(Debug, Clone)]
pub struct ResidueHit<'a> {
    pub comp_id: Option<&'a str>,
    pub chain_id: Option<&'a str>,
    pub seq_id: Option<i32>,
    pub distance: f64,
}

/// Dense distance matrix for a subset of atoms.
#[derive(Debug, Clone)]
pub struct DenseDistanceMatrix {
    /// Atom indices (into `MmcifStructure::atoms`) included in the matrix, in
    /// the same order as the `distances` rows/columns.
    pub atom_indices: Vec<usize>,
    /// Symmetric matrix of raw distances (Å) using minimum-image wrapping if a
    /// unit cell is present.
    pub distances: Vec<Vec<f64>>,
}

/// Neighbor entry for the sparse bond adjacency list.
pub type BondNeighbor = BondEdge;

/// Find the closest atoms around `origin`.
///
/// * `max_distance` – ignore hits farther than this value (Å)
/// * `max_results` – truncate the sorted hit list to this size
pub fn nearest_atoms<'a>(
    structure: &'a MmcifStructure,
    origin: AtomPosition,
    max_distance: f64,
    max_results: usize,
) -> Vec<AtomHit<'a>> {
    if max_distance <= 0.0 || max_results == 0 {
        return Vec::new();
    }
    let unit_cell = structure.metadata.unit_cell;
    let mut hits = Vec::new();
    for (index, atom) in structure.atoms().iter().enumerate() {
        if let Some(pos) = atom.position() {
            let distance = point_distance(origin, pos, unit_cell);
            if distance <= max_distance {
                hits.push(AtomHit {
                    index,
                    distance,
                    atom,
                });
            }
        }
    }
    sort_and_truncate(&mut hits, max_results, |hit| hit.distance);
    hits
}

/// Build a dense distance matrix for the provided `atom_indices`. Atoms without
/// 3D coordinates are skipped. If `atom_indices` is empty, all atoms with known
/// coordinates are used.
pub fn build_dense_distance_matrix(
    structure: &MmcifStructure,
    atom_indices: &[usize],
) -> DenseDistanceMatrix {
    let unit_cell = structure.metadata.unit_cell;
    let mut indices = Vec::new();
    let mut positions = Vec::new();

    if atom_indices.is_empty() {
        for (idx, atom) in structure.atoms().iter().enumerate() {
            if let Some(pos) = atom.position() {
                indices.push(idx);
                positions.push(pos);
            }
        }
    } else {
        for &idx in atom_indices {
            if let Some(atom) = structure.atoms().get(idx) {
                if let Some(pos) = atom.position() {
                    indices.push(idx);
                    positions.push(pos);
                }
            }
        }
    }

    let n = positions.len();
    let mut distances = vec![vec![0.0f64; n]; n];
    match unit_cell {
        Some(cell) => {
            for i in 0..n {
                for j in (i + 1)..n {
                    let distance = cell.minimum_image_distance(positions[i], positions[j]);
                    distances[i][j] = distance;
                    distances[j][i] = distance;
                }
            }
        }
        None => {
            for i in 0..n {
                for j in (i + 1)..n {
                    let distance = positions[i].distance(positions[j]);
                    distances[i][j] = distance;
                    distances[j][i] = distance;
                }
            }
        }
    }

    DenseDistanceMatrix {
        atom_indices: indices,
        distances,
    }
}

/// Find the closest amino-acid residues (based on `_atom_site.label_comp_id`).
pub fn nearest_amino_acids<'a>(
    structure: &'a MmcifStructure,
    origin: AtomPosition,
    max_distance: f64,
    max_results: usize,
) -> Vec<ResidueHit<'a>> {
    if max_distance <= 0.0 || max_results == 0 {
        return Vec::new();
    }
    let unit_cell = structure.metadata.unit_cell;
    let atoms = structure.atoms();
    let residue_keys: Vec<Option<ResidueKey<'a>>> = atoms
        .iter()
        .map(|atom| {
            if is_amino_acid(atom.label_comp_id.as_deref()) {
                Some(ResidueKey::from_atom(atom))
            } else {
                None
            }
        })
        .collect();

    let mut residues: HashMap<ResidueKey<'a>, f64> = HashMap::with_capacity(
        residue_keys
            .iter()
            .filter(|entry| entry.is_some())
            .count()
            .max(8),
    );

    let mut upsert = |key: ResidueKey<'a>, distance: f64| {
        if distance > max_distance {
            return;
        }
        residues
            .entry(key)
            .and_modify(|current| {
                if distance < *current {
                    *current = distance;
                }
            })
            .or_insert(distance);
    };

    match unit_cell {
        Some(cell) => {
            for (atom, key) in atoms.iter().zip(residue_keys.iter()) {
                let Some(key) = key else {
                    continue;
                };
                let Some(pos) = atom.position() else {
                    continue;
                };
                let distance = cell.minimum_image_distance(origin, pos);
                upsert(*key, distance);
            }
        }
        None => {
            for (atom, key) in atoms.iter().zip(residue_keys.iter()) {
                let Some(key) = key else {
                    continue;
                };
                let Some(pos) = atom.position() else {
                    continue;
                };
                let distance = origin.distance(pos);
                upsert(*key, distance);
            }
        }
    }

    let mut hits: Vec<_> = residues
        .into_iter()
        .map(
            |(
                ResidueKey {
                    comp_id,
                    chain_id,
                    seq_id,
                },
                distance,
            )| ResidueHit {
                comp_id,
                chain_id,
                seq_id,
                distance,
            },
        )
        .collect();
    sort_and_truncate(&mut hits, max_results, |hit| hit.distance);
    hits
}

/// Compute a conservative bond list by distance using covalent radii.
///
/// * `tolerance` is an Å value added to the sum of covalent radii to allow
///   for flexibility in real structures (recommended 0.3–0.6 Å).
pub fn compute_bonds_by_distance(
    structure: &MmcifStructure,
    tolerance: f64,
) -> Vec<(usize, usize)> {
    let unit_cell = structure.metadata.unit_cell;
    let atoms = structure.atoms();
    let mut candidates: Vec<(usize, AtomPosition, f64)> = Vec::with_capacity(atoms.len());
    for (index, atom) in atoms.iter().enumerate() {
        if let (Some(position), Some(radius)) = (
            atom.position(),
            atom.element_data().and_then(|e| e.covalent_radius()),
        ) {
            candidates.push((index, position, radius));
        }
    }

    let mut bonds = Vec::new();
    if candidates.len() >= 2 {
        bonds.reserve(candidates.len().saturating_mul(2));
    }

    match unit_cell {
        Some(cell) => {
            for i in 0..candidates.len() {
                let (first_index, first_pos, first_radius) = candidates[i];
                for j in (i + 1)..candidates.len() {
                    let (second_index, second_pos, second_radius) = candidates[j];
                    let distance = cell.minimum_image_distance(first_pos, second_pos);
                    if distance <= first_radius + second_radius + tolerance {
                        bonds.push((first_index, second_index));
                    }
                }
            }
        }
        None => {
            for i in 0..candidates.len() {
                let (first_index, first_pos, first_radius) = candidates[i];
                for j in (i + 1)..candidates.len() {
                    let (second_index, second_pos, second_radius) = candidates[j];
                    let distance = first_pos.distance(second_pos);
                    if distance <= first_radius + second_radius + tolerance {
                        bonds.push((first_index, second_index));
                    }
                }
            }
        }
    }

    bonds
}

/// Resolve atom pairs from explicit bond records already present in the
/// structure.
pub fn resolve_explicit_bonds(structure: &MmcifStructure) -> Vec<ResolvedBond> {
    let mut bonds = Vec::new();
    for (index, neighbors) in structure
        .explicit_bond_adjacency_with_types()
        .iter()
        .enumerate()
    {
        for neighbor in neighbors {
            if index < neighbor.index {
                bonds.push(ResolvedBond {
                    first: index,
                    second: neighbor.index,
                    bond_type: neighbor.bond_type,
                });
            }
        }
    }
    bonds
}

/// Return resolved atom index pairs for explicit mmCIF bond records.
pub fn explicit_bond_pairs(structure: &MmcifStructure) -> Vec<(usize, usize)> {
    resolve_explicit_bonds(structure)
        .into_iter()
        .map(|bond| (bond.first, bond.second))
        .collect()
}

/// Build a sparse adjacency list keyed by bond type. Explicit bonds are used
/// preferentially; optionally, geometric bonds (distance-based) are added as
/// `CovalentUnknown` edges when `include_geometric` is true.
pub fn build_sparse_bond_adjacency(
    structure: &MmcifStructure,
    include_geometric: bool,
    tolerance: f64,
) -> Vec<Vec<BondNeighbor>> {
    let mut adjacency = structure.explicit_bond_adjacency_with_types().to_vec();

    if include_geometric {
        for (i, j) in compute_bonds_by_distance(structure, tolerance) {
            if adjacency[i].iter().any(|edge| edge.index == j) {
                continue;
            }
            adjacency[i].push(BondNeighbor {
                index: j,
                bond_type: Some(BondType::CovalentUnknown),
            });
            adjacency[j].push(BondNeighbor {
                index: i,
                bond_type: Some(BondType::CovalentUnknown),
            });
        }
    }

    adjacency
}

/// Return electronegativity for an atom index if known.
pub fn atom_electronegativity(structure: &MmcifStructure, atom_index: usize) -> Option<f64> {
    structure
        .atoms()
        .get(atom_index)
        .and_then(|a| a.electronegativity())
}

/// Infer a simple hybridization for `atom_index` using a precomputed adjacency
/// list (from `compute_bonds_by_distance`). This uses heavy-atom neighbor
/// counts and is intentionally heuristic.
pub fn infer_hybridization(
    structure: &MmcifStructure,
    atom_index: usize,
    adjacency: &Vec<Vec<usize>>,
) -> Hybridization {
    let atoms = structure.atoms();
    let atom = match atoms.get(atom_index) {
        Some(a) => a,
        None => return Hybridization::Unknown,
    };
    let element = atom.element_data().map(|e| e.symbol).unwrap_or("");
    let neighbors = match adjacency.get(atom_index) {
        Some(list) => list,
        None => return Hybridization::Unknown,
    };
    let total_neighbors = neighbors.len();

    // count heavy neighbors (non-hydrogen)
    let mut heavy_neighbors = 0usize;
    for &nb in neighbors {
        if let Some(ed) = atoms[nb].element_data() {
            if ed.atomic_number > 1 && ed.symbol != "H" {
                heavy_neighbors += 1;
            }
        }
    }

    match element {
        "C" => match total_neighbors {
            4 => Hybridization::Sp3,
            3 => Hybridization::Sp2,
            2 => Hybridization::Sp,
            _ if heavy_neighbors >= 3 => Hybridization::Sp3,
            _ => Hybridization::Unknown,
        },
        "N" => match total_neighbors {
            4 => Hybridization::Sp3,
            3 => Hybridization::Sp2,
            2 => Hybridization::Sp,
            _ if heavy_neighbors >= 3 => Hybridization::Sp3,
            _ => Hybridization::Unknown,
        },
        "O" | "S" => match total_neighbors {
            4 | 3 => Hybridization::Sp3,
            2 => Hybridization::Sp2,
            1 => Hybridization::Sp,
            _ if heavy_neighbors >= 2 => Hybridization::Sp2,
            _ => Hybridization::Unknown,
        },
        _ => Hybridization::Unknown,
    }
}

/// Build an adjacency list from the bond tuple list.
pub fn build_adjacency(natoms: usize, bonds: &[(usize, usize)]) -> Vec<Vec<usize>> {
    if natoms == 0 {
        return Vec::new();
    }

    let mut degrees = vec![0usize; natoms];
    for &(i, j) in bonds {
        if i < natoms && j < natoms {
            degrees[i] += 1;
            degrees[j] += 1;
        }
    }

    let mut adjacency: Vec<Vec<usize>> = degrees
        .into_iter()
        .map(|deg| Vec::with_capacity(deg))
        .collect();

    for &(i, j) in bonds {
        if i < natoms && j < natoms {
            adjacency[i].push(j);
            adjacency[j].push(i);
        }
    }

    adjacency
}

fn point_distance(origin: AtomPosition, target: AtomPosition, unit_cell: Option<UnitCell>) -> f64 {
    if let Some(cell) = unit_cell {
        cell.minimum_image_distance(origin, target)
    } else {
        origin.distance(target)
    }
}

fn sort_and_truncate<T, F>(items: &mut Vec<T>, max_results: usize, mut key: F)
where
    F: FnMut(&T) -> f64,
{
    items.sort_by(|lhs, rhs| match key(lhs).partial_cmp(&key(rhs)) {
        Some(order) => order,
        None => Ordering::Equal,
    });
    if items.len() > max_results {
        items.truncate(max_results);
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct ResidueKey<'a> {
    comp_id: Option<&'a str>,
    chain_id: Option<&'a str>,
    seq_id: Option<i32>,
}

impl<'a> ResidueKey<'a> {
    fn from_atom(atom: &'a AtomSite) -> Self {
        Self {
            comp_id: atom.label_comp_id.as_deref(),
            chain_id: atom
                .label_asym_id
                .as_deref()
                .or(atom.auth_asym_id.as_deref()),
            seq_id: atom.label_seq_id.or(atom.auth_seq_id),
        }
    }
}

fn is_amino_acid(comp_id: Option<&str>) -> bool {
    let Some(raw) = comp_id else {
        return false;
    };
    let upper = raw.trim().to_ascii_uppercase();
    AMINO_ACIDS.iter().any(|&aa| aa == upper)
}

const AMINO_ACIDS: [&str; 21] = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC",
];

#[cfg(test)]
mod tests {
    use super::*;
    use mmcif_core::{
        BondAtomRef, BondRecord, BondRecordSource, BondType, MmcifMetadata, parse_file,
    };
    use std::path::Path;

    fn point(x: f64, y: f64, z: f64) -> AtomPosition {
        AtomPosition { x, y, z }
    }

    fn build_atom(comp_id: &str, chain: &str, seq: i32, coords: (f64, f64, f64)) -> AtomSite {
        let mut atom = AtomSite::new();
        atom.set_value("_atom_site.label_comp_id", comp_id).unwrap();
        atom.set_value("_atom_site.label_asym_id", chain).unwrap();
        atom.set_value("_atom_site.label_seq_id", &seq.to_string())
            .unwrap();
        atom.set_value("_atom_site.Cartn_x", &coords.0.to_string())
            .unwrap();
        atom.set_value("_atom_site.Cartn_y", &coords.1.to_string())
            .unwrap();
        atom.set_value("_atom_site.Cartn_z", &coords.2.to_string())
            .unwrap();
        atom
    }

    fn structure_with_atoms(atoms: Vec<AtomSite>) -> MmcifStructure {
        MmcifStructure::new(MmcifMetadata::default(), atoms, Vec::new())
    }

    #[test]
    fn nearest_atoms_limits_results() {
        let atoms = vec![
            build_atom("GLY", "A", 1, (0.0, 0.0, 0.0)),
            build_atom("GLY", "A", 2, (2.0, 0.0, 0.0)),
            build_atom("GLY", "A", 3, (5.0, 0.0, 0.0)),
        ];
        let structure = structure_with_atoms(atoms);
        let hits = nearest_atoms(&structure, point(0.0, 0.0, 0.0), 3.0, 2);
        assert_eq!(hits.len(), 2);
        assert!(hits[0].distance <= hits[1].distance);
        assert_eq!(hits[0].index, 0);
    }

    #[test]
    fn nearest_residues_filters_amino_acids() {
        let mut atoms = Vec::new();
        atoms.push(build_atom("HOH", "A", 1, (0.0, 0.0, 0.0))); // water, ignored
        atoms.push(build_atom("GLY", "A", 2, (1.0, 0.0, 0.0)));
        atoms.push(build_atom("GLY", "A", 2, (1.5, 0.0, 0.0))); // same residue, farther
        atoms.push(build_atom("ARG", "B", 5, (4.0, 0.0, 0.0)));
        let structure = structure_with_atoms(atoms);
        let hits = nearest_amino_acids(&structure, point(0.0, 0.0, 0.0), 5.0, 5);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].comp_id, Some("GLY"));
        assert!(hits[0].distance < hits[1].distance);
    }

    #[test]
    fn compute_bonds_and_hybridization_basic() {
        // two carbons at ~1.54 Å should be bonded
        let mut a = AtomSite::new();
        a.set_value("_atom_site.type_symbol", "C").unwrap();
        a.set_value("_atom_site.Cartn_x", "0.0").unwrap();
        a.set_value("_atom_site.Cartn_y", "0.0").unwrap();
        a.set_value("_atom_site.Cartn_z", "0.0").unwrap();

        let mut b = AtomSite::new();
        b.set_value("_atom_site.type_symbol", "C").unwrap();
        b.set_value("_atom_site.Cartn_x", "1.54").unwrap();
        b.set_value("_atom_site.Cartn_y", "0.0").unwrap();
        b.set_value("_atom_site.Cartn_z", "0.0").unwrap();

        let structure = MmcifStructure::new(MmcifMetadata::default(), vec![a, b], Vec::new());
        let bonds = compute_bonds_by_distance(&structure, 0.35);
        assert_eq!(bonds.len(), 1);

        // methane-like: central carbon with 4 hydrogens => sp3
        let mut c = AtomSite::new();
        c.set_value("_atom_site.type_symbol", "C").unwrap();
        c.set_value("_atom_site.Cartn_x", "0.0").unwrap();
        c.set_value("_atom_site.Cartn_y", "0.0").unwrap();
        c.set_value("_atom_site.Cartn_z", "0.0").unwrap();

        let mut h1 = AtomSite::new();
        h1.set_value("_atom_site.type_symbol", "H").unwrap();
        h1.set_value("_atom_site.Cartn_x", "0.63").unwrap();
        h1.set_value("_atom_site.Cartn_y", "0.63").unwrap();
        h1.set_value("_atom_site.Cartn_z", "0.63").unwrap();

        let mut h2 = h1.clone();
        h2.set_value("_atom_site.Cartn_x", "-0.63").unwrap();
        let mut h3 = h1.clone();
        h3.set_value("_atom_site.Cartn_y", "-0.63").unwrap();
        let mut h4 = h1.clone();
        h4.set_value("_atom_site.Cartn_z", "-0.63").unwrap();

        let methane = MmcifStructure::new(
            MmcifMetadata::default(),
            vec![c, h1, h2, h3, h4],
            Vec::new(),
        );
        let bonds = compute_bonds_by_distance(&methane, 0.45);
        let adj = build_adjacency(methane.atoms().len(), &bonds);
        let hyb = infer_hybridization(&methane, 0, &adj);
        assert_eq!(hyb, Hybridization::Sp3);
    }

    #[test]
    fn resolve_explicit_bonds_matches_atoms() {
        let mut atom1 = AtomSite::new();
        atom1.set_value("_atom_site.label_atom_id", "CA").unwrap();
        atom1.set_value("_atom_site.label_comp_id", "GLY").unwrap();
        atom1.set_value("_atom_site.label_asym_id", "A").unwrap();
        atom1.set_value("_atom_site.label_seq_id", "1").unwrap();
        atom1.set_value("_atom_site.Cartn_x", "0.0").unwrap();
        atom1.set_value("_atom_site.Cartn_y", "0.0").unwrap();
        atom1.set_value("_atom_site.Cartn_z", "0.0").unwrap();

        let mut atom2 = AtomSite::new();
        atom2.set_value("_atom_site.label_atom_id", "CB").unwrap();
        atom2.set_value("_atom_site.label_comp_id", "GLY").unwrap();
        atom2.set_value("_atom_site.label_asym_id", "A").unwrap();
        atom2.set_value("_atom_site.label_seq_id", "1").unwrap();
        atom2.set_value("_atom_site.Cartn_x", "1.5").unwrap();
        atom2.set_value("_atom_site.Cartn_y", "0.0").unwrap();
        atom2.set_value("_atom_site.Cartn_z", "0.0").unwrap();

        let mut structure = MmcifStructure::default();
        structure.atoms.push(atom1);
        structure.atoms.push(atom2);

        let mut record = BondRecord::new(BondRecordSource::StructConn);
        record.atom1.label_atom_id = Some("CA".to_string());
        record.atom1.label_comp_id = Some("GLY".to_string());
        record.atom1.label_asym_id = Some("A".to_string());
        record.atom1.label_seq_id = Some(1);
        record.atom2.label_atom_id = Some("CB".to_string());
        record.atom2.label_comp_id = Some("GLY".to_string());
        record.atom2.label_asym_id = Some("A".to_string());
        record.atom2.label_seq_id = Some(1);
        record.order = Some("SING".to_string());
        record.update_bond_type_hint();
        structure.bonds.push(record);

        let resolved = resolve_explicit_bonds(&structure);
        assert_eq!(resolved.len(), 1);
        assert_eq!(resolved[0].first, 0);
        assert_eq!(resolved[0].second, 1);
        assert_eq!(resolved[0].bond_type, Some(BondType::CovalentSingle));
        let pair_list: Vec<(usize, usize)> = resolved.iter().map(|b| (b.first, b.second)).collect();
        let adjacency = build_adjacency(structure.atoms().len(), &pair_list);
        assert_eq!(adjacency[0], vec![1]);
    }

    #[test]
    fn explicit_bond_pairs_map_to_indices() {
        let mut a = build_atom("ALA", "A", 1, (0.0, 0.0, 0.0));
        a.set_value("_atom_site.label_atom_id", "CA").unwrap();
        let mut b = build_atom("ALA", "A", 1, (1.5, 0.0, 0.0));
        b.set_value("_atom_site.label_atom_id", "N").unwrap();

        let mut atom1 = BondAtomRef::default();
        atom1.label_atom_id = Some("CA".to_string());
        atom1.label_comp_id = Some("ALA".to_string());
        atom1.label_asym_id = Some("A".to_string());
        atom1.label_seq_id = Some(1);
        let mut atom2 = BondAtomRef::default();
        atom2.label_atom_id = Some("N".to_string());
        atom2.label_comp_id = Some("ALA".to_string());
        atom2.label_asym_id = Some("A".to_string());
        atom2.label_seq_id = Some(1);
        let mut bond = BondRecord::new(BondRecordSource::StructConn);
        bond.atom1 = atom1;
        bond.atom2 = atom2;
        bond.connection_type = Some("covalent".to_string());
        bond.order = Some("DOUB".to_string());
        bond.distance = Some(1.5);
        bond.update_bond_type_hint();

        let structure = MmcifStructure::new(MmcifMetadata::default(), vec![a, b], vec![bond]);

        let resolved = resolve_explicit_bonds(&structure);
        assert_eq!(resolved[0].bond_type, Some(BondType::CovalentDouble));
        assert_eq!(explicit_bond_pairs(&structure), vec![(0, 1)]);
    }

    #[test]
    fn dense_distance_matrix_subset() {
        let atoms = vec![
            build_atom("GLY", "A", 1, (0.0, 0.0, 0.0)),
            build_atom("GLY", "A", 1, (3.0, 0.0, 0.0)),
            build_atom("GLY", "A", 1, (0.0, 4.0, 0.0)),
        ];
        let structure = structure_with_atoms(atoms);
        let matrix = build_dense_distance_matrix(&structure, &[0, 2]);
        assert_eq!(matrix.atom_indices, vec![0, 2]);
        assert_eq!(matrix.distances.len(), 2);
        assert!((matrix.distances[0][1] - 4.0).abs() < 1e-6);
        assert_eq!(matrix.distances[0][0], 0.0);
    }

    #[test]
    fn sparse_bond_adjacency_combines_sources() {
        let mut atoms = vec![
            build_atom("GLY", "A", 1, (0.0, 0.0, 0.0)),
            build_atom("GLY", "A", 1, (1.4, 0.0, 0.0)),
            build_atom("GLY", "A", 1, (2.8, 0.0, 0.0)),
        ];
        for atom in &mut atoms {
            atom.set_value("_atom_site.type_symbol", "C").unwrap();
        }

        let mut record = BondRecord::new(BondRecordSource::StructConn);
        record.atom1.label_atom_id = Some("CA".to_string());
        record.atom1.label_comp_id = Some("GLY".to_string());
        record.atom1.label_asym_id = Some("A".to_string());
        record.atom1.label_seq_id = Some(1);
        record.atom2.label_atom_id = Some("N".to_string());
        record.atom2.label_comp_id = Some("GLY".to_string());
        record.atom2.label_asym_id = Some("A".to_string());
        record.atom2.label_seq_id = Some(1);
        record.order = Some("DOUB".to_string());
        record.update_bond_type_hint();

        let mut structure = MmcifStructure::new(MmcifMetadata::default(), atoms, vec![record]);

        // Ensure atom labels align with references
        structure.atoms[0]
            .set_value("_atom_site.label_atom_id", "CA")
            .unwrap();
        structure.atoms[0]
            .set_value("_atom_site.label_comp_id", "GLY")
            .unwrap();
        structure.atoms[0]
            .set_value("_atom_site.label_asym_id", "A")
            .unwrap();
        structure.atoms[0]
            .set_value("_atom_site.label_seq_id", "1")
            .unwrap();
        structure.atoms[1]
            .set_value("_atom_site.label_atom_id", "N")
            .unwrap();
        structure.atoms[1]
            .set_value("_atom_site.label_comp_id", "GLY")
            .unwrap();
        structure.atoms[1]
            .set_value("_atom_site.label_asym_id", "A")
            .unwrap();
        structure.atoms[1]
            .set_value("_atom_site.label_seq_id", "1")
            .unwrap();

        let adjacency = build_sparse_bond_adjacency(&structure, true, 0.35);
        assert_eq!(adjacency.len(), 3);
        let first_neighbors = &adjacency[0];
        assert_eq!(first_neighbors.len(), 1);
        assert_eq!(
            first_neighbors[0],
            BondNeighbor {
                index: 1,
                bond_type: Some(BondType::CovalentDouble),
            }
        );

        let second_neighbors = &adjacency[1];
        assert!(
            second_neighbors
                .iter()
                .any(|n| n.index == 0 && n.bond_type == Some(BondType::CovalentDouble))
        );
        assert!(
            second_neighbors
                .iter()
                .any(|n| n.index == 2 && n.bond_type == Some(BondType::CovalentUnknown))
        );
    }

    #[test]
    fn integration_real_structure_11as() {
        let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
        let cif_path = manifest.join("../../11AS.cif.gz");
        let structure = parse_file(&cif_path).expect("parse 11AS structure");
        assert!(structure.atom_count() > 0);

        let n_idx = structure
            .atoms()
            .iter()
            .position(|atom| matches_atom(atom, "ALA", "A", "N"))
            .expect("ALA A N present");
        let anchor_atom = &structure.atoms()[n_idx];
        let origin = anchor_atom.position().expect("N atom coordinates");
        let seq_id = anchor_atom
            .label_seq_id
            .or(anchor_atom.auth_seq_id)
            .expect("sequence id on anchor atom");
        let chain_id = anchor_atom
            .label_asym_id
            .as_deref()
            .or(anchor_atom.auth_asym_id.as_deref())
            .expect("chain id on anchor atom");
        let comp_id = anchor_atom
            .label_comp_id
            .as_deref()
            .expect("component id on anchor atom");
        let ca_idx = find_atom_in_residue(&structure, comp_id, chain_id, seq_id, "CA");
        let c_idx = find_atom_in_residue(&structure, comp_id, chain_id, seq_id, "C");

        let hits = nearest_atoms(&structure, origin, 2.5, 6);
        assert_eq!(hits[0].index, n_idx);
        assert!(hits[0].distance <= f64::EPSILON);
        let ca_hit = hits
            .iter()
            .find(|hit| hit.index == ca_idx)
            .expect("CA near N");
        assert_close(ca_hit.distance, 1.4914922058, 1e-6);
        let c_hit = hits
            .iter()
            .find(|hit| hit.index == c_idx)
            .expect("C near N");
        assert_close(c_hit.distance, 2.4604505685, 1e-6);

        let residue_hits = nearest_amino_acids(&structure, origin, 3.0, 1);
        assert!(!residue_hits.is_empty());
        assert_eq!(residue_hits[0].comp_id, Some("ALA"));
        assert_eq!(residue_hits[0].chain_id, Some("A"));
        assert_eq!(residue_hits[0].seq_id, Some(seq_id));

        let explicit = resolve_explicit_bonds(&structure);
        assert!(
            explicit
                .iter()
                .any(|bond| matches_pair(bond, n_idx, ca_idx))
        );

        let adjacency = build_sparse_bond_adjacency(&structure, true, 0.35);
        let n_neighbors = &adjacency[n_idx];
        assert!(n_neighbors.iter().any(|neighbor| {
            neighbor.index == ca_idx && neighbor.bond_type == Some(BondType::CovalentSingle)
        }));

        let matrix = build_dense_distance_matrix(&structure, &[n_idx, ca_idx, c_idx]);
        assert_eq!(matrix.atom_indices, vec![n_idx, ca_idx, c_idx]);
        assert_close(matrix.distances[0][1], 1.4914922058126883, 1e-6);
        assert_close(matrix.distances[0][2], 2.460450568493501, 1e-6);
        assert_close(matrix.distances[1][2], 1.53347513837036, 1e-6);
    }

    fn find_atom_in_residue(
        structure: &MmcifStructure,
        comp_id: &str,
        chain_id: &str,
        seq_id: i32,
        atom_id: &str,
    ) -> usize {
        structure
            .atoms()
            .iter()
            .position(|atom| {
                atom.label_comp_id.as_deref() == Some(comp_id)
                    && atom
                        .label_asym_id
                        .as_deref()
                        .or(atom.auth_asym_id.as_deref())
                        == Some(chain_id)
                    && atom.label_seq_id.or(atom.auth_seq_id) == Some(seq_id)
                    && atom.label_atom_id.as_deref() == Some(atom_id)
            })
            .expect("atom present in structure")
    }

    fn matches_atom(atom: &AtomSite, comp_id: &str, chain_id: &str, atom_id: &str) -> bool {
        atom.label_comp_id.as_deref() == Some(comp_id)
            && atom
                .label_asym_id
                .as_deref()
                .or(atom.auth_asym_id.as_deref())
                == Some(chain_id)
            && atom.label_atom_id.as_deref() == Some(atom_id)
    }

    fn matches_pair(bond: &ResolvedBond, a: usize, b: usize) -> bool {
        (bond.first == a && bond.second == b) || (bond.first == b && bond.second == a)
    }

    fn assert_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "value mismatch: got {}, expected {} (tol {})",
            actual,
            expected,
            tolerance
        );
    }
}
