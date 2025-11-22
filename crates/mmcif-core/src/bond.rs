use crate::atom::{AtomSite, BondType};

/// A reference to an atom in the mmCIF data, as provided by `_struct_conn` or
/// `_chem_comp_bond` records.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BondAtomRef {
    pub label_atom_id: Option<String>,
    pub label_comp_id: Option<String>,
    pub label_asym_id: Option<String>,
    pub label_seq_id: Option<i32>,
    pub auth_atom_id: Option<String>,
    pub auth_asym_id: Option<String>,
    pub auth_seq_id: Option<i32>,
}

impl BondAtomRef {
    pub fn matches(&self, atom: &AtomSite) -> bool {
        fn matches_str(expected: &Option<String>, actual: Option<&str>) -> bool {
            match expected {
                Some(expected) => actual
                    .map(|value| value.eq_ignore_ascii_case(expected))
                    .unwrap_or(false),
                None => true,
            }
        }

        fn matches_seq(expected: Option<i32>, actual: Option<i32>) -> bool {
            match expected {
                Some(expected) => actual.map_or(false, |value| value == expected),
                None => true,
            }
        }

        matches_str(&self.label_atom_id, atom.label_atom_id.as_deref())
            && matches_str(&self.auth_atom_id, atom.auth_atom_id.as_deref())
            && matches_str(&self.label_comp_id, atom.label_comp_id.as_deref())
            && matches_str(&self.label_asym_id, atom.label_asym_id.as_deref())
            && matches_str(&self.auth_asym_id, atom.auth_asym_id.as_deref())
            && matches_seq(self.label_seq_id, atom.label_seq_id)
            && matches_seq(self.auth_seq_id, atom.auth_seq_id)
    }
}

/// Source tag for parsed bond records.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondRecordSource {
    StructConn,
    ChemCompBond,
    Unknown,
}

impl Default for BondRecordSource {
    fn default() -> Self {
        BondRecordSource::Unknown
    }
}

/// Parsed representation of an explicit bond record from mmCIF.
#[derive(Clone, Debug, PartialEq)]
pub struct BondRecord {
    pub source: BondRecordSource,
    pub atom1: BondAtomRef,
    pub atom2: BondAtomRef,
    pub connection_type: Option<String>,
    pub order: Option<String>,
    pub distance: Option<f64>,
    pub bond_type: Option<BondType>,
}

impl BondRecord {
    pub fn new(source: BondRecordSource) -> Self {
        Self {
            source,
            atom1: BondAtomRef::default(),
            atom2: BondAtomRef::default(),
            connection_type: None,
            order: None,
            distance: None,
            bond_type: None,
        }
    }

    pub fn update_bond_type_hint(&mut self) {
        if let Some(order) = self.order.as_deref() {
            if let Some(kind) = BondType::from_order_str(order) {
                self.bond_type = Some(kind);
                return;
            }
        }

        if let Some(conn) = self.connection_type.as_deref() {
            if let Some(kind) = BondType::from_connection_type(conn) {
                self.bond_type = Some(kind);
            }
        }
    }

    pub fn bond_type(&self) -> Option<BondType> {
        self.bond_type
    }
}

impl Default for BondRecord {
    fn default() -> Self {
        Self::new(BondRecordSource::Unknown)
    }
}

/// Edge entry in an adjacency list referencing a bonded neighbor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BondEdge {
    pub index: usize,
    pub bond_type: Option<BondType>,
}
