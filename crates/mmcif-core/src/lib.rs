#![forbid(unsafe_code)]
//! Core primitives for reading mmCIF/mmCIF PDBx data files.
//! The modules intentionally mirror the structure of `pdbtbx` and
//! Biopython's `MMCIFParser`, but the implementation is local so this
//! workspace can evolve independently.

pub mod atom;
pub mod bond;
pub mod error;
pub mod model;
pub mod parser;

pub use atom::{AtomPosition, AtomSite, BondType, ContactClass, ElementData};
pub use bond::{BondAtomRef, BondEdge, BondRecord, BondRecordSource};
pub use error::{ParseError, ParseErrorKind};
pub use model::{InterAtomicContact, MmcifMetadata, MmcifStructure, UnitCell};
pub use parser::{ParserOptions, parse_async_reader, parse_file, parse_reader, parse_str};
