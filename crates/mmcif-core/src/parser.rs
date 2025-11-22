use std::{
    fs::File,
    io::{BufRead, Read},
    path::Path,
};

use flate2::read::GzDecoder;
use futures::io::{AsyncRead, AsyncReadExt};

use crate::{
    atom::AtomSite,
    bond::{BondRecord, BondRecordSource},
    error::{ParseError, ParseErrorKind},
    model::{MmcifMetadata, MmcifStructure, UnitCell},
};

#[derive(Clone, Copy, Debug)]
pub struct ParserOptions {
    pub require_atom_site_loop: bool,
}

impl Default for ParserOptions {
    fn default() -> Self {
        Self {
            require_atom_site_loop: true,
        }
    }
}

pub fn parse_file(path: impl AsRef<Path>) -> Result<MmcifStructure, ParseError> {
    let path_ref = path.as_ref();
    let file = File::open(path_ref).map_err(|err| ParseError::from(err).with_path(path_ref))?;
    let mut reader: Box<dyn Read> = if path_ref
        .extension()
        .and_then(|ext| ext.to_str())
        .map(|ext| ext.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut buf = String::new();
    reader.read_to_string(&mut buf)?;
    parse_str_with_options(&buf, ParserOptions::default()).map(|mut structure| {
        if structure.metadata.entry_id.is_none() {
            structure.metadata.entry_id = path_ref
                .file_stem()
                .and_then(|stem| stem.to_str())
                .map(|s| s.to_string());
        }
        structure
    })
}

pub fn parse_reader<R: BufRead>(
    mut reader: R,
    options: ParserOptions,
) -> Result<MmcifStructure, ParseError> {
    let mut buf = String::new();
    reader.read_to_string(&mut buf)?;
    parse_str_with_options(&buf, options)
}

pub async fn parse_async_reader<R>(
    mut reader: R,
    options: ParserOptions,
) -> Result<MmcifStructure, ParseError>
where
    R: AsyncRead + Unpin,
{
    let mut buf = Vec::new();
    reader
        .read_to_end(&mut buf)
        .await
        .map_err(ParseError::from)?;
    let text = String::from_utf8(buf)
        .map_err(|_| ParseError::new(ParseErrorKind::Utf8, "input is not valid UTF-8"))?;
    parse_str_with_options(&text, options)
}

pub fn parse_str(input: &str) -> Result<MmcifStructure, ParseError> {
    parse_str_with_options(input, ParserOptions::default())
}

fn parse_str_with_options(
    input: &str,
    options: ParserOptions,
) -> Result<MmcifStructure, ParseError> {
    let tokens = tokenize(input)?;
    Parser::new(tokens, options).run()
}

struct Parser {
    tokens: Vec<String>,
    cursor: usize,
    options: ParserOptions,
    metadata: MmcifMetadata,
    atoms: Vec<AtomSite>,
    unit_cell_builder: UnitCellBuilder,
    bonds: Vec<BondRecord>,
}

impl Parser {
    fn new(tokens: Vec<String>, options: ParserOptions) -> Self {
        Self {
            tokens,
            cursor: 0,
            options,
            metadata: MmcifMetadata::default(),
            atoms: Vec::new(),
            unit_cell_builder: UnitCellBuilder::default(),
            bonds: Vec::new(),
        }
    }

    fn run(mut self) -> Result<MmcifStructure, ParseError> {
        while let Some(token) = self.peek().cloned() {
            if token.eq_ignore_ascii_case("loop_") {
                self.cursor += 1;
                self.parse_loop()?;
            } else if token.starts_with("data_") {
                self.cursor += 1;
                self.metadata.entry_id = Some(token[5..].to_string());
            } else if token.starts_with('_') {
                self.cursor += 1;
                let value = self.next_value()?;
                self.assign_scalar(&token, &value)?;
            } else {
                self.cursor += 1; // Ignore other tokens for now
            }
        }

        if self.metadata.unit_cell.is_none() {
            self.metadata.unit_cell = self.unit_cell_builder.build()?;
        }

        if self.options.require_atom_site_loop && self.atoms.is_empty() {
            return Err(ParseError::new(
                ParseErrorKind::MissingAtomLoop,
                "_atom_site loop missing",
            ));
        }

        Ok(MmcifStructure::new(self.metadata, self.atoms, self.bonds))
    }

    fn parse_loop(&mut self) -> Result<(), ParseError> {
        let tags = self.collect_loop_tags()?;
        if tags.is_empty() {
            return Err(ParseError::new(
                ParseErrorKind::LoopMismatch,
                "loop_ without tags",
            ));
        }

        let is_atom_loop = tags
            .iter()
            .all(|tag| tag.to_ascii_lowercase().starts_with("_atom_site."));
        let is_struct_conn_loop = tags
            .iter()
            .any(|tag| tag.to_ascii_lowercase().starts_with("_struct_conn."));
        let is_chem_comp_loop = tags
            .iter()
            .any(|tag| tag.to_ascii_lowercase().starts_with("_chem_comp_bond."));

        let mut values = Vec::new();
        while let Some(token) = self.peek() {
            if token.eq_ignore_ascii_case("loop_")
                || token.starts_with('_')
                || token.starts_with("data_")
            {
                break;
            }
            values.push(token.clone());
            self.cursor += 1;
        }

        if values.len() % tags.len() != 0 {
            if is_atom_loop {
                return Err(ParseError::new(
                    ParseErrorKind::LoopMismatch,
                    format!(
                        "loop contains {} values but {} tags",
                        values.len(),
                        tags.len()
                    ),
                ));
            }
            return Ok(());
        }

        if is_atom_loop {
            for chunk in values.chunks(tags.len()) {
                let mut site = AtomSite::new();
                for (tag, value) in tags.iter().zip(chunk.iter()) {
                    site.set_value(tag, value)?;
                }
                self.atoms.push(site);
            }
        } else if is_struct_conn_loop {
            self.parse_struct_conn_loop(&tags, &values)?;
        } else if is_chem_comp_loop {
            self.parse_chem_comp_bond_loop(&tags, &values)?;
        }

        Ok(())
    }

    fn parse_struct_conn_loop(
        &mut self,
        tags: &[String],
        values: &[String],
    ) -> Result<(), ParseError> {
        for chunk in values.chunks(tags.len()) {
            let mut record = BondRecord::new(BondRecordSource::StructConn);
            for (tag, value) in tags.iter().zip(chunk.iter()) {
                let key = tag.to_ascii_lowercase();
                match key.as_str() {
                    "_struct_conn.ptnr1_label_atom_id" => {
                        record.atom1.label_atom_id = normalize_scalar(value)
                    }
                    "_struct_conn.ptnr1_auth_atom_id" => {
                        record.atom1.auth_atom_id = normalize_scalar(value)
                    }
                    "_struct_conn.ptnr1_label_comp_id" => {
                        record.atom1.label_comp_id = normalize_scalar(value)
                    }
                    "_struct_conn.ptnr1_label_asym_id" => {
                        record.atom1.label_asym_id = normalize_scalar(value)
                    }
                    "_struct_conn.ptnr1_auth_asym_id" => {
                        record.atom1.auth_asym_id = normalize_scalar(value)
                    }
                    "_struct_conn.ptnr1_label_seq_id" => {
                        record.atom1.label_seq_id = parse_optional_int(value, tag)?;
                    }
                    "_struct_conn.ptnr1_auth_seq_id" => {
                        record.atom1.auth_seq_id = parse_optional_int(value, tag)?;
                    }
                    "_struct_conn.ptnr2_label_atom_id" => {
                        record.atom2.label_atom_id = normalize_scalar(value)
                    }
                    "_struct_conn.ptnr2_auth_atom_id" => {
                        record.atom2.auth_atom_id = normalize_scalar(value)
                    }
                    "_struct_conn.ptnr2_label_comp_id" => {
                        record.atom2.label_comp_id = normalize_scalar(value)
                    }
                    "_struct_conn.ptnr2_label_asym_id" => {
                        record.atom2.label_asym_id = normalize_scalar(value)
                    }
                    "_struct_conn.ptnr2_auth_asym_id" => {
                        record.atom2.auth_asym_id = normalize_scalar(value)
                    }
                    "_struct_conn.ptnr2_label_seq_id" => {
                        record.atom2.label_seq_id = parse_optional_int(value, tag)?;
                    }
                    "_struct_conn.ptnr2_auth_seq_id" => {
                        record.atom2.auth_seq_id = parse_optional_int(value, tag)?;
                    }
                    "_struct_conn.conn_type_id" => {
                        record.connection_type = normalize_scalar(value);
                        record.update_bond_type_hint();
                    }
                    "_struct_conn.pdbx_dist_value" => {
                        record.distance = parse_optional_float(value, tag)?;
                    }
                    _ => {}
                }
            }
            self.bonds.push(record);
        }
        Ok(())
    }

    fn parse_chem_comp_bond_loop(
        &mut self,
        tags: &[String],
        values: &[String],
    ) -> Result<(), ParseError> {
        for chunk in values.chunks(tags.len()) {
            let mut record = BondRecord::new(BondRecordSource::ChemCompBond);
            for (tag, value) in tags.iter().zip(chunk.iter()) {
                let key = tag.to_ascii_lowercase();
                match key.as_str() {
                    "_chem_comp_bond.atom_id_1" => {
                        record.atom1.label_atom_id = normalize_scalar(value)
                    }
                    "_chem_comp_bond.atom_id_2" => {
                        record.atom2.label_atom_id = normalize_scalar(value)
                    }
                    "_chem_comp_bond.comp_id" => {
                        let comp = normalize_scalar(value);
                        record.atom1.label_comp_id = comp.clone();
                        record.atom2.label_comp_id = comp;
                    }
                    "_chem_comp_bond.comp_id_1" => {
                        record.atom1.label_comp_id = normalize_scalar(value)
                    }
                    "_chem_comp_bond.comp_id_2" => {
                        record.atom2.label_comp_id = normalize_scalar(value)
                    }
                    "_chem_comp_bond.value_order" => {
                        record.order = normalize_scalar(value);
                        record.update_bond_type_hint();
                    }
                    _ => {}
                }
            }
            record.update_bond_type_hint();
            self.bonds.push(record);
        }
        Ok(())
    }

    fn collect_loop_tags(&mut self) -> Result<Vec<String>, ParseError> {
        let mut tags = Vec::new();
        while let Some(token) = self.peek() {
            if token.starts_with('_') {
                tags.push(token.clone());
                self.cursor += 1;
            } else {
                break;
            }
        }
        Ok(tags)
    }

    fn next_value(&mut self) -> Result<String, ParseError> {
        match self.peek().cloned() {
            Some(value) => {
                self.cursor += 1;
                Ok(value)
            }
            None => Err(ParseError::new(
                ParseErrorKind::MissingValue,
                "expected value, found EOF",
            )),
        }
    }

    fn peek(&self) -> Option<&String> {
        self.tokens.get(self.cursor)
    }

    fn assign_scalar(&mut self, tag: &str, raw_value: &str) -> Result<(), ParseError> {
        match tag {
            "_entry.id" => self.metadata.entry_id = normalize_scalar(raw_value),
            "_struct.title" => self.metadata.title = normalize_scalar(raw_value),
            "_symmetry.space_group_name_H-M" => {
                self.metadata.space_group = normalize_scalar(raw_value)
            }
            tag if tag.starts_with("_cell.") => self.unit_cell_builder.set(tag, raw_value)?,
            _ => {}
        }
        Ok(())
    }
}

#[derive(Default)]
struct UnitCellBuilder {
    a: Option<f64>,
    b: Option<f64>,
    c: Option<f64>,
    alpha: Option<f64>,
    beta: Option<f64>,
    gamma: Option<f64>,
}

impl UnitCellBuilder {
    fn set(&mut self, tag: &str, value: &str) -> Result<(), ParseError> {
        match tag {
            "_cell.length_a" => self.a = Some(parse_float(value, tag)?),
            "_cell.length_b" => self.b = Some(parse_float(value, tag)?),
            "_cell.length_c" => self.c = Some(parse_float(value, tag)?),
            "_cell.angle_alpha" => self.alpha = Some(parse_float(value, tag)?),
            "_cell.angle_beta" => self.beta = Some(parse_float(value, tag)?),
            "_cell.angle_gamma" => self.gamma = Some(parse_float(value, tag)?),
            _ => {}
        }
        Ok(())
    }

    fn build(self) -> Result<Option<UnitCell>, ParseError> {
        let UnitCellBuilder {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        } = self;
        let mut missing = 0;
        if a.is_none() {
            missing += 1;
        }
        if b.is_none() {
            missing += 1;
        }
        if c.is_none() {
            missing += 1;
        }
        if alpha.is_none() {
            missing += 1;
        }
        if beta.is_none() {
            missing += 1;
        }
        if gamma.is_none() {
            missing += 1;
        }

        if missing == 0 {
            Ok(Some(UnitCell::new(
                a.unwrap(),
                b.unwrap(),
                c.unwrap(),
                alpha.unwrap(),
                beta.unwrap(),
                gamma.unwrap(),
            )?))
        } else if missing == 6 {
            Ok(None)
        } else {
            Err(ParseError::new(
                ParseErrorKind::MissingValue,
                "partial unit cell definition encountered",
            ))
        }
    }
}

fn parse_float(value: &str, tag: &str) -> Result<f64, ParseError> {
    let trimmed = value.trim();
    if trimmed == "?" || trimmed == "." || trimmed.is_empty() {
        return Err(ParseError::new(
            ParseErrorKind::MissingValue,
            format!("required numeric value missing for {tag}"),
        ));
    }
    let cleaned = strip_esd(trimmed);
    cleaned.parse::<f64>().map_err(|_| {
        ParseError::new(
            ParseErrorKind::InvalidNumber,
            format!("invalid float '{cleaned}' for {tag}"),
        )
    })
}

fn parse_optional_int(value: &str, tag: &str) -> Result<Option<i32>, ParseError> {
    if let Some(clean) = normalize_scalar(value) {
        clean.parse::<i32>().map(Some).map_err(|_| {
            ParseError::new(
                ParseErrorKind::InvalidNumber,
                format!("invalid integer '{clean}' for {tag}"),
            )
        })
    } else {
        Ok(None)
    }
}

fn parse_optional_float(value: &str, tag: &str) -> Result<Option<f64>, ParseError> {
    if let Some(clean) = normalize_scalar(value) {
        let stripped = strip_esd(&clean);
        stripped.parse::<f64>().map(Some).map_err(|_| {
            ParseError::new(
                ParseErrorKind::InvalidNumber,
                format!("invalid numeric literal '{stripped}' for {tag}"),
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

fn tokenize(input: &str) -> Result<Vec<String>, ParseError> {
    let mut tokens = Vec::new();
    let bytes = input.as_bytes();
    let mut index = 0;
    let mut line_start = true;

    while index < bytes.len() {
        match bytes[index] {
            b'\n' => {
                line_start = true;
                index += 1;
            }
            b'\r' | b'\t' | b' ' => {
                index += 1;
            }
            b'#' => {
                while index < bytes.len() && bytes[index] != b'\n' {
                    index += 1;
                }
            }
            b';' if line_start => {
                let start = index + 1;
                index += 1;
                let mut end = None;
                while index < bytes.len() {
                    if bytes[index] == b'\n' {
                        if index + 1 < bytes.len() && bytes[index + 1] == b';' {
                            end = Some(index);
                            index += 2; // skip newline + semicolon
                            while index < bytes.len() && bytes[index] != b'\n' {
                                index += 1;
                            }
                            break;
                        }
                    }
                    index += 1;
                }
                let end = end.ok_or_else(|| {
                    ParseError::new(ParseErrorKind::InvalidToken, "unterminated text field")
                })?;
                let text = input[start..end].trim_end_matches(['\r', '\n']).to_string();
                tokens.push(text);
                line_start = true;
            }
            b'"' | b'\'' => {
                let quote = bytes[index];
                index += 1;
                let start = index;
                while index < bytes.len() && bytes[index] != quote {
                    index += 1;
                }
                if index >= bytes.len() {
                    return Err(ParseError::new(
                        ParseErrorKind::InvalidToken,
                        "unterminated quoted value",
                    ));
                }
                tokens.push(input[start..index].to_string());
                index += 1;
                line_start = false;
            }
            _ => {
                let start = index;
                while index < bytes.len() && !matches!(bytes[index], b' ' | b'\t' | b'\r' | b'\n') {
                    if bytes[index] == b'#' {
                        break;
                    }
                    index += 1;
                }
                tokens.push(input[start..index].to_string());
                line_start = false;
            }
        }
    }

    Ok(tokens)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::BondType;
    use crate::bond::BondRecordSource;
    use pretty_assertions::assert_eq;
    use std::path::Path;

    #[test]
    fn parse_minimal_structure() {
        let cif = r#"
        data_demo
        _entry.id demo
        _cell.length_a 10.0
        _cell.length_b 11.0
        _cell.length_c 12.0
        _cell.angle_alpha 90
        _cell.angle_beta 90
        _cell.angle_gamma 90
        loop_
        _atom_site.id
        _atom_site.type_symbol
        _atom_site.Cartn_x
        _atom_site.Cartn_y
        _atom_site.Cartn_z
        1 C 0.0 0.0 0.0
        2 O 1.0 0.0 0.0
        "#;

        let structure = parse_str(cif).expect("parsed");
        assert_eq!(structure.atom_count(), 2);
        assert_eq!(structure.metadata.entry_id, Some("demo".to_string()));
        assert!(structure.metadata.unit_cell.is_some());
    }

    #[test]
    fn parse_real_cif_gz() {
        let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
        let cif_path = manifest.join("../../11AS.cif.gz");
        let structure = parse_file(&cif_path).expect("parse compressed cif");
        assert!(structure.atom_count() > 0);
        assert!(structure.metadata.entry_id.is_some());
    }

    #[test]
    fn parse_struct_conn_bonds() {
        let cif = r#"
        data_bonds
        loop_
        _atom_site.id
        _atom_site.label_atom_id
        _atom_site.label_comp_id
        _atom_site.label_asym_id
        _atom_site.label_seq_id
        1 CA ALA A 1
        2 N ALA A 1
        loop_
        _struct_conn.ptnr1_label_atom_id
        _struct_conn.ptnr1_label_comp_id
        _struct_conn.ptnr1_label_asym_id
        _struct_conn.ptnr1_label_seq_id
        _struct_conn.ptnr2_label_atom_id
        _struct_conn.ptnr2_label_comp_id
        _struct_conn.ptnr2_label_asym_id
        _struct_conn.ptnr2_label_seq_id
        _struct_conn.conn_type_id
        CA ALA A 1 N ALA A 1 covalent
        "#;

        let structure = parse_str(cif).expect("bond loop");
        assert_eq!(structure.bond_records().len(), 1);
        let record = &structure.bond_records()[0];
        assert_eq!(record.source, BondRecordSource::StructConn);
        assert_eq!(record.connection_type.as_deref(), Some("covalent"));
        assert_eq!(structure.resolve_bond(record), Some((0, 1)));
    }

    #[test]
    fn parse_chem_comp_bond_orders() {
        let cif = r#"
        data_comp
        loop_
        _atom_site.id
        _atom_site.label_atom_id
        _atom_site.label_comp_id
        _atom_site.label_asym_id
        _atom_site.label_seq_id
        1 C1 BENZ A 1
        2 C2 BENZ A 1
        loop_
        _chem_comp_bond.atom_id_1
        _chem_comp_bond.atom_id_2
        _chem_comp_bond.comp_id
        _chem_comp_bond.value_order
        C1 C2 BENZ DOUB
        "#;

        let structure = parse_str(cif).expect("chem comp bond");
        assert_eq!(structure.bond_records().len(), 1);
        let record = &structure.bond_records()[0];
        assert_eq!(record.order.as_deref(), Some("DOUB"));
        assert_eq!(record.bond_type(), Some(BondType::CovalentDouble));
    }
}
