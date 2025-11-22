use std::cmp::Ordering;

use crate::error::{ParseError, ParseErrorKind};
use crate::model::UnitCell;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AtomPosition {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl AtomPosition {
    pub fn distance(self, other: Self) -> f64 {
        self.squared_distance(other).sqrt()
    }

    pub fn squared_distance(self, other: Self) -> f64 {
        (other.z - self.z).mul_add(
            other.z - self.z,
            (other.y - self.y).mul_add(other.y - self.y, (other.x - self.x).powi(2)),
        )
    }

    pub fn vector_to(self, other: Self) -> [f64; 3] {
        [other.x - self.x, other.y - self.y, other.z - self.z]
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AtomSite {
    pub id: Option<String>,
    pub label_atom_id: Option<String>,
    pub auth_atom_id: Option<String>,
    pub label_comp_id: Option<String>,
    pub alt_id: Option<String>,
    pub type_symbol: Option<String>,
    pub label_asym_id: Option<String>,
    pub auth_asym_id: Option<String>,
    pub label_seq_id: Option<i32>,
    pub auth_seq_id: Option<i32>,
    pub model_number: Option<i32>,
    pub occupancy: Option<f64>,
    pub b_iso: Option<f64>,
    pub formal_charge: Option<i32>,
    pub cartn_x: Option<f64>,
    pub cartn_y: Option<f64>,
    pub cartn_z: Option<f64>,
}

impl AtomSite {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_value(&mut self, tag: &str, raw_value: &str) -> Result<(), ParseError> {
        match tag {
            "_atom_site.id" => self.id = normalize_scalar(raw_value),
            "_atom_site.label_atom_id" => self.label_atom_id = normalize_scalar(raw_value),
            "_atom_site.auth_atom_id" => self.auth_atom_id = normalize_scalar(raw_value),
            "_atom_site.label_comp_id" => self.label_comp_id = normalize_scalar(raw_value),
            "_atom_site.alt_id" | "_atom_site.label_alt_id" => {
                self.alt_id = normalize_scalar(raw_value)
            }
            "_atom_site.type_symbol" => {
                self.type_symbol = normalize_scalar(raw_value).map(|s| s.to_ascii_uppercase())
            }
            "_atom_site.label_asym_id" => self.label_asym_id = normalize_scalar(raw_value),
            "_atom_site.auth_asym_id" => self.auth_asym_id = normalize_scalar(raw_value),
            "_atom_site.label_seq_id" => self.label_seq_id = parse_optional_int(raw_value, tag)?,
            "_atom_site.auth_seq_id" => self.auth_seq_id = parse_optional_int(raw_value, tag)?,
            "_atom_site.pdbx_PDB_model_num" => {
                self.model_number = parse_optional_int(raw_value, tag)?
            }
            "_atom_site.Cartn_x" => self.cartn_x = parse_optional_float(raw_value, tag)?,
            "_atom_site.Cartn_y" => self.cartn_y = parse_optional_float(raw_value, tag)?,
            "_atom_site.Cartn_z" => self.cartn_z = parse_optional_float(raw_value, tag)?,
            "_atom_site.occupancy" => self.occupancy = parse_optional_float(raw_value, tag)?,
            "_atom_site.B_iso_or_equiv" => self.b_iso = parse_optional_float(raw_value, tag)?,
            "_atom_site.pdbx_formal_charge" => {
                self.formal_charge = parse_optional_int(raw_value, tag)?
            }
            _ => {}
        }
        Ok(())
    }

    pub fn position(&self) -> Option<AtomPosition> {
        Some(AtomPosition {
            x: self.cartn_x?,
            y: self.cartn_y?,
            z: self.cartn_z?,
        })
    }

    pub fn distance_to(&self, other: &Self) -> Option<f64> {
        Some(self.position()?.distance(other.position()?))
    }

    pub fn wrapped_distance_to(&self, other: &Self, cell: &UnitCell) -> Option<f64> {
        let a = self.position()?;
        let b = other.position()?;
        Some(cell.minimum_image_distance(a, b))
    }

    pub fn element_data(&self) -> Option<&'static ElementData> {
        self.type_symbol
            .as_deref()
            .and_then(ElementData::from_symbol)
            .or_else(|| {
                self.label_atom_id
                    .as_deref()
                    .and_then(ElementData::from_symbol)
            })
            .or_else(|| {
                self.auth_atom_id
                    .as_deref()
                    .and_then(ElementData::from_symbol)
            })
    }

    /// Return Pauling electronegativity for this atom if known.
    pub fn electronegativity(&self) -> Option<f64> {
        self.element_data().and_then(|e| e.electronegativity())
    }

    /// Infer a coarse bond type between `self` and `other` using element
    /// information and the observed `distance` (Å). This is heuristic and
    /// conservative – prefer explicit bond records when available.
    pub fn infer_bond_type(&self, other: &Self, distance: f64) -> BondType {
        // prefer covalent classification from radii
        let class = self.classify_contact(other, distance);
        if class == ContactClass::Covalent {
            // try to separate single vs double by comparing to covalent radii sum
            if let (Some(a), Some(b)) = (
                self.element_data().and_then(|e| e.covalent_radius()),
                other.element_data().and_then(|e| e.covalent_radius()),
            ) {
                let expected = a + b;
                // heuristics: significantly shorter than expected -> tighter/multiple bond
                if distance < expected - 0.15 {
                    return BondType::CovalentDouble;
                }
                return BondType::CovalentSingle;
            }
            return BondType::CovalentUnknown;
        }

        // hydrogen-bond heuristic
        if let (Some(ea), Some(eb)) = (self.element_data(), other.element_data()) {
            let a_sym = ea.symbol;
            let b_sym = eb.symbol;
            if (a_sym == "H" && (b_sym == "O" || b_sym == "N"))
                || (b_sym == "H" && (a_sym == "O" || a_sym == "N"))
            {
                if distance > 1.2 && distance < 2.8 {
                    return BondType::HydrogenBond;
                }
            }
            // ionic heuristic: large electronegativity difference
            if let (Some(en_a), Some(en_b)) = (ea.electronegativity, eb.electronegativity) {
                let diff = (en_a - en_b).abs();
                if diff >= 1.7 {
                    return BondType::Ionic;
                }
            }
        }

        BondType::Unknown
    }

    pub fn classify_contact(&self, other: &Self, distance: f64) -> ContactClass {
        let element_a = self.element_data();
        let element_b = other.element_data();
        classify_by_radii(distance, element_a, element_b)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ContactClass {
    Covalent,
    ShortContact,
    Clash,
    NonBonded,
}

impl ContactClass {
    pub fn severity(&self) -> Ordering {
        match self {
            ContactClass::Covalent => Ordering::Equal,
            ContactClass::ShortContact => Ordering::Less,
            ContactClass::Clash => Ordering::Greater,
            ContactClass::NonBonded => Ordering::Equal,
        }
    }
}

/// Rough bond type classification returned by inference heuristics.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondType {
    CovalentSingle,
    CovalentDouble,
    CovalentTriple,
    CovalentAromatic,
    CovalentUnknown,
    Ionic,
    HydrogenBond,
    Metallic,
    Unknown,
}

impl BondType {
    pub fn from_order_str(order: &str) -> Option<Self> {
        let upper = order.trim().to_ascii_uppercase();
        match upper.as_str() {
            "1" | "SING" | "SINGLE" => Some(BondType::CovalentSingle),
            "2" | "DOUB" | "DOUBLE" => Some(BondType::CovalentDouble),
            "3" | "TRIP" | "TRIPLE" => Some(BondType::CovalentTriple),
            "4" | "QUAD" | "QUADRUPLE" => Some(BondType::CovalentUnknown),
            "AROM" | "AROMATIC" | "DELO" | "DELOCALIZED" => {
                Some(BondType::CovalentAromatic)
            }
            "PI" => Some(BondType::CovalentDouble),
            "PARTIAL" => Some(BondType::CovalentUnknown),
            "ION" | "IONIC" => Some(BondType::Ionic),
            "METAL" | "COORD" => Some(BondType::Metallic),
            _ => None,
        }
    }

    pub fn from_connection_type(conn: &str) -> Option<Self> {
        let upper = conn.trim().to_ascii_uppercase();
        if upper.contains("HYD") {
            return Some(BondType::HydrogenBond);
        }
        if upper.contains("ION") {
            return Some(BondType::Ionic);
        }
        if upper.contains("METAL") {
            return Some(BondType::Metallic);
        }
        if upper.contains("COVAL") {
            return Some(BondType::CovalentUnknown);
        }
        None
    }
}

/// Simple hybridization guesses for common organic atoms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Hybridization {
    Sp,
    Sp2,
    Sp3,
    Aromatic,
    Unknown,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElementData {
    pub symbol: &'static str,
    pub atomic_number: u8,
    pub covalent_radius: Option<f64>,
    pub van_der_waals_radius: Option<f64>,
    pub electronegativity: Option<f64>,
}

impl ElementData {
    pub const fn new(
        symbol: &'static str,
        atomic_number: u8,
        covalent_radius: Option<f64>,
        van_der_waals_radius: Option<f64>,
        electronegativity: Option<f64>,
    ) -> Self {
        Self {
            symbol,
            atomic_number,
            covalent_radius,
            van_der_waals_radius,
            electronegativity,
        }
    }

    pub fn from_symbol(symbol: &str) -> Option<&'static Self> {
        let upper = symbol.trim().to_ascii_uppercase();
        ELEMENT_TABLE
            .iter()
            .find(|element| element.symbol.eq_ignore_ascii_case(&upper))
    }

    /// Return Pauling electronegativity if available.
    pub fn electronegativity(&self) -> Option<f64> {
        self.electronegativity
    }

    /// Return covalent radius in Å if available.
    pub fn covalent_radius(&self) -> Option<f64> {
        self.covalent_radius
    }
}

const ELEMENT_TABLE: [ElementData; 18] = [
    ElementData::new("H", 1, Some(0.32), Some(1.2), Some(2.20)),
    ElementData::new("C", 6, Some(0.75), Some(1.7), Some(2.55)),
    ElementData::new("N", 7, Some(0.71), Some(1.55), Some(3.04)),
    ElementData::new("O", 8, Some(0.63), Some(1.52), Some(3.44)),
    ElementData::new("F", 9, Some(0.64), Some(1.47), Some(3.98)),
    ElementData::new("P", 15, Some(1.06), Some(1.8), Some(2.19)),
    ElementData::new("S", 16, Some(1.02), Some(1.8), Some(2.58)),
    ElementData::new("CL", 17, Some(0.99), Some(1.75), Some(3.16)),
    ElementData::new("NA", 11, Some(1.66), Some(2.27), Some(0.93)),
    ElementData::new("MG", 12, Some(1.41), Some(1.73), Some(1.31)),
    ElementData::new("K", 19, Some(2.03), Some(2.75), Some(0.82)),
    ElementData::new("CA", 20, Some(1.74), Some(2.31), Some(1.00)),
    ElementData::new("FE", 26, Some(1.16), Some(1.63), Some(1.83)),
    ElementData::new("ZN", 30, Some(1.25), Some(1.39), Some(1.65)),
    ElementData::new("SE", 34, Some(1.17), Some(1.9), Some(2.55)),
    ElementData::new("BR", 35, Some(1.21), Some(1.85), Some(2.96)),
    ElementData::new("I", 53, Some(1.40), Some(1.98), Some(2.66)),
    ElementData::new("CU", 29, Some(1.22), Some(1.4), Some(1.90)),
];

fn classify_by_radii(
    distance: f64,
    lhs: Option<&ElementData>,
    rhs: Option<&ElementData>,
) -> ContactClass {
    let target = lhs.and_then(|l| rhs.map(|r| (l, r))).and_then(|(l, r)| {
        match (l.covalent_radius, r.covalent_radius) {
            (Some(a), Some(b)) => Some(a + b),
            _ => None,
        }
    });

    if let Some(expected) = target {
        if (distance - expected).abs() <= 0.4 {
            return ContactClass::Covalent;
        }
        if distance < expected - 0.25 {
            return ContactClass::Clash;
        }
        if distance < expected + 1.0 {
            return ContactClass::ShortContact;
        }
    } else if distance < 2.0 {
        return ContactClass::ShortContact;
    }

    ContactClass::NonBonded
}

fn parse_optional_float(value: &str, tag: &str) -> Result<Option<f64>, ParseError> {
    if let Some(clean) = normalize_scalar(value) {
        let trimmed = strip_esd(&clean);
        trimmed.parse::<f64>().map(Some).map_err(|_| {
            ParseError::new(
                ParseErrorKind::InvalidNumber,
                format!("invalid numeric literal '{trimmed}' for {tag}"),
            )
        })
    } else {
        Ok(None)
    }
}

fn parse_optional_int(value: &str, tag: &str) -> Result<Option<i32>, ParseError> {
    if let Some(clean) = normalize_scalar(value) {
        clean.parse::<i32>().map(Some).map_err(|_| {
            ParseError::new(
                ParseErrorKind::InvalidNumber,
                format!("invalid integer literal '{clean}' for {tag}"),
            )
        })
    } else {
        Ok(None)
    }
}

fn normalize_scalar(value: &str) -> Option<String> {
    let trimmed = value.trim();
    if trimmed.is_empty() || trimmed == "?" || trimmed == "." {
        None
    } else {
        Some(trimmed.to_string())
    }
}

fn strip_esd(value: &str) -> &str {
    if let Some(idx) = value.find('(') {
        &value[..idx]
    } else {
        value
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn parse_atom_site_values() {
        let mut site = AtomSite::new();
        site.set_value("_atom_site.id", "1").unwrap();
        site.set_value("_atom_site.type_symbol", "c").unwrap();
        site.set_value("_atom_site.Cartn_x", "10.0").unwrap();
        site.set_value("_atom_site.Cartn_y", "11.0").unwrap();
        site.set_value("_atom_site.Cartn_z", "12.0").unwrap();
        assert_eq!(site.type_symbol.as_deref(), Some("C"));
        assert!(site.position().is_some());
    }

    #[test]
    fn classify_distance() {
        let mut a = AtomSite::new();
        let mut b = AtomSite::new();
        a.set_value("_atom_site.type_symbol", "C").unwrap();
        b.set_value("_atom_site.type_symbol", "C").unwrap();
        assert_eq!(a.classify_contact(&b, 1.54), ContactClass::Covalent);
        assert_eq!(a.classify_contact(&b, 3.5), ContactClass::NonBonded);
    }
}
