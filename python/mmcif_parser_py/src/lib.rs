#![forbid(unsafe_code)]

use std::path::PathBuf;

use mmcif_analyze::{
    atom_electronegativity, build_dense_distance_matrix, build_sparse_bond_adjacency,
    explicit_bond_pairs, infer_hybridization, nearest_amino_acids, nearest_atoms, BondNeighbor,
};
use mmcif_core::atom::Hybridization;
use mmcif_core::{parse_file, AtomPosition, BondType, ContactClass, MmcifStructure};
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyAny;
use thiserror::Error;

#[derive(Error, Debug)]
enum BindingError {
    #[error("failed to read mmCIF: {0}")]
    Parse(#[from] mmcif_core::ParseError),
    #[error("path must be a string or os.PathLike")]
    Path,
}

impl From<BindingError> for PyErr {
    fn from(value: BindingError) -> Self {
        match value {
            BindingError::Parse(err) => PyValueError::new_err(err.to_string()),
            BindingError::Path => PyIOError::new_err("path must be str or os.PathLike"),
        }
    }
}

#[pyclass(module = "mmcif_parser", name = "Structure")]
pub struct PyStructure {
    inner: MmcifStructure,
}

#[pymethods]
impl PyStructure {
    #[getter]
    pub fn atom_count(&self) -> usize {
        self.inner.atom_count()
    }

    #[getter]
    pub fn entry_id(&self) -> Option<String> {
        self.inner.metadata.entry_id.clone()
    }

    #[getter]
    pub fn space_group(&self) -> Option<String> {
        self.inner.metadata.space_group.clone()
    }

    pub fn unit_cell(&self) -> Option<((f64, f64, f64), (f64, f64, f64))> {
        self.inner.metadata.unit_cell.map(|cell| {
            let lengths = cell.lengths();
            let angles = cell.angles();
            (lengths, angles)
        })
    }

    pub fn bounding_box(&self) -> Option<((f64, f64, f64), (f64, f64, f64))> {
        self.inner
            .bounding_box()
            .map(|(min, max)| ((min.x, min.y, min.z), (max.x, max.y, max.z)))
    }

    pub fn contacts(&self, cutoff: f64, max_results: Option<usize>) -> PyResult<Vec<PyContact>> {
        if cutoff <= 0.0 {
            return Err(PyValueError::new_err("cutoff must be positive"));
        }
        let mut contacts: Vec<PyContact> = self
            .inner
            .inspect_contacts(cutoff)
            .into_iter()
            .map(PyContact::from)
            .collect();
        contacts.sort_by(|lhs, rhs| {
            lhs.distance
                .partial_cmp(&rhs.distance)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        if let Some(limit) = max_results {
            if contacts.len() > limit {
                contacts.truncate(limit);
            }
        }
        Ok(contacts)
    }

    pub fn nearest_atoms(
        &self,
        origin: (f64, f64, f64),
        max_distance: f64,
        max_results: usize,
    ) -> Vec<PyAtomHit> {
        let hits = nearest_atoms(&self.inner, to_position(origin), max_distance, max_results);
        hits.into_iter().map(PyAtomHit::from).collect()
    }

    pub fn nearest_residues(
        &self,
        origin: (f64, f64, f64),
        max_distance: f64,
        max_results: usize,
    ) -> Vec<PyResidueHit> {
        let hits = nearest_amino_acids(&self.inner, to_position(origin), max_distance, max_results);
        hits.into_iter().map(PyResidueHit::from).collect()
    }

    pub fn dense_distance_matrix(&self, atom_indices: Option<Vec<usize>>) -> PyDistanceMatrix {
        let indices = atom_indices.unwrap_or_default();
        let matrix = build_dense_distance_matrix(&self.inner, &indices);
        PyDistanceMatrix {
            atom_indices: matrix.atom_indices,
            distances: matrix.distances,
        }
    }

    pub fn sparse_adjacency(
        &self,
        include_geometric: bool,
        tolerance: f64,
    ) -> Vec<Vec<PyBondNeighbor>> {
        build_sparse_bond_adjacency(&self.inner, include_geometric, tolerance)
            .into_iter()
            .map(|neighbors| neighbors.into_iter().map(PyBondNeighbor::from).collect())
            .collect()
    }

    pub fn explicit_bonds(&self) -> Vec<(usize, usize)> {
        explicit_bond_pairs(&self.inner)
    }

    pub fn atom_profile(
        &self,
        atom_index: usize,
        include_geometric: bool,
        tolerance: f64,
    ) -> PyResult<PyAtomProfile> {
        if atom_index >= self.inner.atom_count() {
            return Err(PyValueError::new_err("atom_index out of range"));
        }
        if include_geometric && tolerance < 0.0 {
            return Err(PyValueError::new_err("tolerance must be non-negative"));
        }
        let typed_adjacency: Vec<Vec<BondNeighbor>> = if include_geometric {
            build_sparse_bond_adjacency(&self.inner, true, tolerance)
        } else {
            self.inner
                .explicit_bond_adjacency_with_types()
                .iter()
                .cloned()
                .collect()
        };
        let adjacency: Vec<Vec<usize>> = typed_adjacency
            .iter()
            .map(|neighbors| neighbors.iter().map(|edge| edge.index).collect())
            .collect();

        let hybridization = infer_hybridization(&self.inner, atom_index, &adjacency);
        let electronegativity = atom_electronegativity(&self.inner, atom_index);
        let atom = &self.inner.atoms()[atom_index];
        let chain_id = atom
            .label_asym_id
            .clone()
            .or_else(|| atom.auth_asym_id.clone());
        let seq_id = atom.label_seq_id.or(atom.auth_seq_id);
        let position = atom.position().map(|pos| (pos.x, pos.y, pos.z));

        Ok(PyAtomProfile {
            index: atom_index,
            label_atom_id: atom.label_atom_id.clone(),
            label_comp_id: atom.label_comp_id.clone(),
            chain_id,
            seq_id,
            element: atom.type_symbol.clone(),
            occupancy: atom.occupancy,
            b_factor: atom.b_iso,
            formal_charge: atom.formal_charge,
            position,
            electronegativity,
            hybridization: hybridization_label(hybridization).map(str::to_string),
            neighbor_indices: adjacency.get(atom_index).cloned().unwrap_or_default(),
            neighbor_bonds: typed_adjacency
                .get(atom_index)
                .cloned()
                .unwrap_or_default()
                .into_iter()
                .map(PyBondNeighbor::from)
                .collect(),
        })
    }
}

#[pyclass(module = "mmcif_parser", name = "AtomHit")]
#[derive(Clone)]
pub struct PyAtomHit {
    #[pyo3(get)]
    pub index: usize,
    #[pyo3(get)]
    pub distance: f64,
    #[pyo3(get)]
    pub label_atom_id: Option<String>,
    #[pyo3(get)]
    pub label_comp_id: Option<String>,
    #[pyo3(get)]
    pub chain_id: Option<String>,
    #[pyo3(get)]
    pub seq_id: Option<i32>,
    #[pyo3(get)]
    pub element: Option<String>,
    #[pyo3(get)]
    pub position: Option<(f64, f64, f64)>,
}

impl From<mmcif_analyze::AtomHit<'_>> for PyAtomHit {
    fn from(value: mmcif_analyze::AtomHit<'_>) -> Self {
        let atom = value.atom;
        let position = atom.position().map(|pos| (pos.x, pos.y, pos.z));
        Self {
            index: value.index,
            distance: value.distance,
            label_atom_id: atom.label_atom_id.clone(),
            label_comp_id: atom.label_comp_id.clone(),
            chain_id: atom
                .label_asym_id
                .clone()
                .or_else(|| atom.auth_asym_id.clone()),
            seq_id: atom.label_seq_id.or(atom.auth_seq_id),
            element: atom.type_symbol.clone(),
            position,
        }
    }
}

#[pyclass(module = "mmcif_parser", name = "ResidueHit")]
#[derive(Clone)]
pub struct PyResidueHit {
    #[pyo3(get)]
    pub comp_id: Option<String>,
    #[pyo3(get)]
    pub chain_id: Option<String>,
    #[pyo3(get)]
    pub seq_id: Option<i32>,
    #[pyo3(get)]
    pub distance: f64,
}

impl From<mmcif_analyze::ResidueHit<'_>> for PyResidueHit {
    fn from(value: mmcif_analyze::ResidueHit<'_>) -> Self {
        Self {
            comp_id: value.comp_id.map(|s| s.to_string()),
            chain_id: value.chain_id.map(|s| s.to_string()),
            seq_id: value.seq_id,
            distance: value.distance,
        }
    }
}

#[pyclass(module = "mmcif_parser", name = "DistanceMatrix")]
pub struct PyDistanceMatrix {
    #[pyo3(get)]
    pub atom_indices: Vec<usize>,
    #[pyo3(get)]
    pub distances: Vec<Vec<f64>>,
}

#[pyclass(module = "mmcif_parser", name = "Contact")]
pub struct PyContact {
    #[pyo3(get)]
    pub first: usize,
    #[pyo3(get)]
    pub second: usize,
    #[pyo3(get)]
    pub distance: f64,
    #[pyo3(get)]
    pub classification: String,
}

impl From<mmcif_core::InterAtomicContact> for PyContact {
    fn from(value: mmcif_core::InterAtomicContact) -> Self {
        Self {
            first: value.first,
            second: value.second,
            distance: value.distance,
            classification: contact_class_label(value.class).to_string(),
        }
    }
}

#[pyclass(module = "mmcif_parser", name = "BondNeighbor")]
#[derive(Clone)]
pub struct PyBondNeighbor {
    #[pyo3(get)]
    pub index: usize,
    #[pyo3(get)]
    pub bond_type: Option<String>,
}

impl From<BondNeighbor> for PyBondNeighbor {
    fn from(value: BondNeighbor) -> Self {
        Self {
            index: value.index,
            bond_type: value.bond_type.map(bond_type_label),
        }
    }
}

#[pyclass(module = "mmcif_parser", name = "AtomProfile")]
pub struct PyAtomProfile {
    #[pyo3(get)]
    pub index: usize,
    #[pyo3(get)]
    pub label_atom_id: Option<String>,
    #[pyo3(get)]
    pub label_comp_id: Option<String>,
    #[pyo3(get)]
    pub chain_id: Option<String>,
    #[pyo3(get)]
    pub seq_id: Option<i32>,
    #[pyo3(get)]
    pub element: Option<String>,
    #[pyo3(get)]
    pub occupancy: Option<f64>,
    #[pyo3(get)]
    pub b_factor: Option<f64>,
    #[pyo3(get)]
    pub formal_charge: Option<i32>,
    #[pyo3(get)]
    pub position: Option<(f64, f64, f64)>,
    #[pyo3(get)]
    pub electronegativity: Option<f64>,
    #[pyo3(get)]
    pub hybridization: Option<String>,
    #[pyo3(get)]
    pub neighbor_indices: Vec<usize>,
    #[pyo3(get)]
    pub neighbor_bonds: Vec<PyBondNeighbor>,
}

fn bond_type_label(bond: BondType) -> String {
    match bond {
        BondType::CovalentSingle => "covalent_single",
        BondType::CovalentDouble => "covalent_double",
        BondType::CovalentTriple => "covalent_triple",
        BondType::CovalentAromatic => "covalent_aromatic",
        BondType::CovalentUnknown => "covalent_unknown",
        BondType::HydrogenBond => "hydrogen_bond",
        BondType::Ionic => "ionic",
        BondType::Metallic => "metallic",
        BondType::Unknown => "unknown",
    }
    .to_string()
}

fn contact_class_label(class: ContactClass) -> &'static str {
    match class {
        ContactClass::Covalent => "covalent",
        ContactClass::ShortContact => "short_contact",
        ContactClass::Clash => "clash",
        ContactClass::NonBonded => "nonbonded",
    }
}

fn hybridization_label(value: Hybridization) -> Option<&'static str> {
    match value {
        Hybridization::Sp => Some("sp"),
        Hybridization::Sp2 => Some("sp2"),
        Hybridization::Sp3 => Some("sp3"),
        Hybridization::Aromatic => Some("aromatic"),
        Hybridization::Unknown => None,
    }
}

#[pyfunction]
pub fn parse(path: &Bound<'_, PyAny>) -> PyResult<PyStructure> {
    let path_buf: PathBuf = path.extract().map_err(|_| BindingError::Path)?;
    let structure = parse_file(&path_buf).map_err(BindingError::from)?;
    Ok(PyStructure { inner: structure })
}

#[pymodule]
fn mmcif_parser(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyStructure>()?;
    m.add_class::<PyAtomHit>()?;
    m.add_class::<PyResidueHit>()?;
    m.add_class::<PyDistanceMatrix>()?;
    m.add_class::<PyBondNeighbor>()?;
    m.add_class::<PyContact>()?;
    m.add_class::<PyAtomProfile>()?;
    m.add_function(wrap_pyfunction!(parse, m)?)?;
    Ok(())
}

fn to_position(origin: (f64, f64, f64)) -> AtomPosition {
    AtomPosition {
        x: origin.0,
        y: origin.1,
        z: origin.2,
    }
}
