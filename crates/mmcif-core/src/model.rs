use crate::atom::{AtomPosition, AtomSite, ContactClass};
use crate::bond::{BondAtomRef, BondEdge, BondRecord};
use crate::error::{ParseError, ParseErrorKind};
use once_cell::unsync::OnceCell;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct MmcifMetadata {
    pub entry_id: Option<String>,
    pub title: Option<String>,
    pub space_group: Option<String>,
    pub unit_cell: Option<UnitCell>,
}

#[derive(Debug, Default)]
pub struct MmcifStructure {
    pub metadata: MmcifMetadata,
    pub atoms: Vec<AtomSite>,
    pub bonds: Vec<BondRecord>,
    cache: BondAdjacencyCache,
}

impl MmcifStructure {
    pub fn new(metadata: MmcifMetadata, atoms: Vec<AtomSite>, bonds: Vec<BondRecord>) -> Self {
        Self {
            metadata,
            atoms,
            bonds,
            cache: BondAdjacencyCache::default(),
        }
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn atoms(&self) -> &[AtomSite] {
        &self.atoms
    }

    pub fn atoms_mut(&mut self) -> &mut [AtomSite] {
        self.clear_caches();
        &mut self.atoms
    }

    pub fn bond_records(&self) -> &[BondRecord] {
        &self.bonds
    }

    pub fn bonds_mut(&mut self) -> &mut [BondRecord] {
        self.clear_caches();
        &mut self.bonds
    }

    pub fn clear_caches(&mut self) {
        self.cache.clear();
    }

    pub fn find_atom_index(&self, reference: &BondAtomRef) -> Option<usize> {
        self.atoms.iter().position(|atom| reference.matches(atom))
    }

    pub fn resolve_bond(&self, record: &BondRecord) -> Option<(usize, usize)> {
        let first = self.find_atom_index(&record.atom1)?;
        let second = self.find_atom_index(&record.atom2)?;
        Some((first, second))
    }

    pub fn bounding_box(&self) -> Option<(AtomPosition, AtomPosition)> {
        let mut iter = self.atoms.iter().filter_map(AtomSite::position);
        let first = iter.next()?;
        let mut min = first;
        let mut max = first;
        for pos in iter {
            min.x = min.x.min(pos.x);
            min.y = min.y.min(pos.y);
            min.z = min.z.min(pos.z);
            max.x = max.x.max(pos.x);
            max.y = max.y.max(pos.y);
            max.z = max.z.max(pos.z);
        }
        Some((min, max))
    }

    pub fn inspect_contacts(&self, cutoff: f64) -> Vec<InterAtomicContact> {
        assert!(cutoff > 0.0, "cutoff must be positive");
        let mut contacts = Vec::new();
        let unit_cell = self.metadata.unit_cell;
        for (i, atom_a) in self.atoms.iter().enumerate() {
            let pos_a = match atom_a.position() {
                Some(pos) => pos,
                None => continue,
            };
            for (j, atom_b) in self.atoms.iter().enumerate().skip(i + 1) {
                let pos_b = match atom_b.position() {
                    Some(pos) => pos,
                    None => continue,
                };
                let mut distance = pos_a.distance(pos_b);
                if let Some(cell) = unit_cell {
                    distance = cell.minimum_image_distance(pos_a, pos_b);
                }
                if distance <= cutoff {
                    let class = atom_a.classify_contact(atom_b, distance);
                    contacts.push(InterAtomicContact {
                        first: i,
                        second: j,
                        distance,
                        class,
                    });
                }
            }
        }
        contacts
    }

    pub fn explicit_bond_adjacency(&self) -> &[Vec<usize>] {
        self.cache
            .untyped
            .get_or_init(|| self.build_explicit_bond_adjacency())
    }

    pub fn explicit_bond_adjacency_with_types(&self) -> &[Vec<BondEdge>] {
        self.cache
            .typed
            .get_or_init(|| self.build_typed_explicit_bond_adjacency())
    }

    fn build_explicit_bond_adjacency(&self) -> Vec<Vec<usize>> {
        let mut adjacency = vec![Vec::new(); self.atoms.len()];
        for record in &self.bonds {
            if let Some((first, second)) = self.resolve_bond(record) {
                adjacency[first].push(second);
                adjacency[second].push(first);
            }
        }
        adjacency
    }

    fn build_typed_explicit_bond_adjacency(&self) -> Vec<Vec<BondEdge>> {
        let mut adjacency = vec![Vec::new(); self.atoms.len()];
        for record in &self.bonds {
            if let Some((first, second)) = self.resolve_bond(record) {
                adjacency[first].push(BondEdge {
                    index: second,
                    bond_type: record.bond_type(),
                });
                adjacency[second].push(BondEdge {
                    index: first,
                    bond_type: record.bond_type(),
                });
            }
        }
        adjacency
    }
}

impl Clone for MmcifStructure {
    fn clone(&self) -> Self {
        Self {
            metadata: self.metadata.clone(),
            atoms: self.atoms.clone(),
            bonds: self.bonds.clone(),
            cache: BondAdjacencyCache::default(),
        }
    }
}

impl PartialEq for MmcifStructure {
    fn eq(&self, other: &Self) -> bool {
        self.metadata == other.metadata && self.atoms == other.atoms && self.bonds == other.bonds
    }
}

#[derive(Default, Debug)]
struct BondAdjacencyCache {
    untyped: OnceCell<Vec<Vec<usize>>>,
    typed: OnceCell<Vec<Vec<BondEdge>>>,
}

impl BondAdjacencyCache {
    fn clear(&mut self) {
        self.untyped = OnceCell::new();
        self.typed = OnceCell::new();
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UnitCell {
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
}

impl UnitCell {
    pub fn new(
        a: f64,
        b: f64,
        c: f64,
        alpha: f64,
        beta: f64,
        gamma: f64,
    ) -> Result<Self, ParseError> {
        if !(a > 0.0 && b > 0.0 && c > 0.0) {
            return Err(ParseError::new(
                ParseErrorKind::InvalidNumber,
                "unit cell lengths must be positive",
            ));
        }
        if !(valid_angle(alpha) && valid_angle(beta) && valid_angle(gamma)) {
            return Err(ParseError::new(
                ParseErrorKind::InvalidNumber,
                "unit cell angles must be within (0, 180)",
            ));
        }
        Ok(Self {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        })
    }

    #[inline]
    pub fn lengths(&self) -> (f64, f64, f64) {
        (self.a, self.b, self.c)
    }

    #[inline]
    pub fn angles(&self) -> (f64, f64, f64) {
        (self.alpha, self.beta, self.gamma)
    }

    #[inline]
    pub fn a(&self) -> f64 {
        self.a
    }

    #[inline]
    pub fn b(&self) -> f64 {
        self.b
    }

    #[inline]
    pub fn c(&self) -> f64 {
        self.c
    }

    #[inline]
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    #[inline]
    pub fn beta(&self) -> f64 {
        self.beta
    }

    #[inline]
    pub fn gamma(&self) -> f64 {
        self.gamma
    }

    pub fn minimum_image_distance(&self, first: AtomPosition, second: AtomPosition) -> f64 {
        if self.is_orthogonal() {
            let dx = wrap_delta(second.x - first.x, self.a);
            let dy = wrap_delta(second.y - first.y, self.b);
            let dz = wrap_delta(second.z - first.z, self.c);
            AtomPosition {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            }
            .distance(AtomPosition {
                x: dx,
                y: dy,
                z: dz,
            })
        } else {
            first.distance(second)
        }
    }

    fn is_orthogonal(&self) -> bool {
        (self.alpha - 90.0).abs() < 1e-4
            && (self.beta - 90.0).abs() < 1e-4
            && (self.gamma - 90.0).abs() < 1e-4
    }
}

fn wrap_delta(delta: f64, length: f64) -> f64 {
    let half = length / 2.0;
    if delta > half {
        delta - length
    } else if delta < -half {
        delta + length
    } else {
        delta
    }
}

fn valid_angle(angle: f64) -> bool {
    angle > 0.0 && angle < 180.0
}

#[derive(Clone, Debug, PartialEq)]
pub struct InterAtomicContact {
    pub first: usize,
    pub second: usize,
    pub distance: f64,
    pub class: ContactClass,
}

impl InterAtomicContact {
    pub fn atoms<'a>(&self, structure: &'a MmcifStructure) -> Option<(&'a AtomSite, &'a AtomSite)> {
        let a = structure.atoms.get(self.first)?;
        let b = structure.atoms.get(self.second)?;
        Some((a, b))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::AtomSite;
    use pretty_assertions::assert_eq;

    #[test]
    fn bounding_box() {
        let mut structure = MmcifStructure::default();
        let mut atom = AtomSite::new();
        atom.set_value("_atom_site.Cartn_x", "0.0").unwrap();
        atom.set_value("_atom_site.Cartn_y", "1.0").unwrap();
        atom.set_value("_atom_site.Cartn_z", "2.0").unwrap();
        structure.atoms.push(atom.clone());
        atom.set_value("_atom_site.Cartn_x", "2.0").unwrap();
        atom.set_value("_atom_site.Cartn_y", "3.0").unwrap();
        atom.set_value("_atom_site.Cartn_z", "4.0").unwrap();
        structure.atoms.push(atom);
        let bbox = structure.bounding_box().unwrap();
        assert!((bbox.0.x - 0.0).abs() < 1e-6);
        assert!((bbox.1.z - 4.0).abs() < 1e-6);
    }

    #[test]
    fn detect_contacts() {
        let mut structure = MmcifStructure::default();
        let mut atom_a = AtomSite::new();
        atom_a.set_value("_atom_site.type_symbol", "C").unwrap();
        atom_a.set_value("_atom_site.Cartn_x", "0.0").unwrap();
        atom_a.set_value("_atom_site.Cartn_y", "0.0").unwrap();
        atom_a.set_value("_atom_site.Cartn_z", "0.0").unwrap();
        let mut atom_b = AtomSite::new();
        atom_b.set_value("_atom_site.type_symbol", "C").unwrap();
        atom_b.set_value("_atom_site.Cartn_x", "1.5").unwrap();
        atom_b.set_value("_atom_site.Cartn_y", "0.0").unwrap();
        atom_b.set_value("_atom_site.Cartn_z", "0.0").unwrap();
        structure.atoms.push(atom_a);
        structure.atoms.push(atom_b);
        let contacts = structure.inspect_contacts(2.0);
        assert_eq!(contacts.len(), 1);
        assert_eq!(contacts[0].class, ContactClass::Covalent);
    }
}
