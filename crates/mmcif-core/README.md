# mmcif-core

Lightweight primitives to parse mmCIF/mmCIF PDBx files and expose strongly typed atom records. The code mirrors the project layout of `pdbtbx`/Biopython, but every module is implemented locally so the workspace can evolve independently and stay dependency-free.

## Implementation highlights

- **Tokeniser + tolerant loop parser**: `parser::tokenize` performs the minimal lexical split required for mmCIF (quoted values, multiline `;` records, comments). During parsing we only enforce tag/value arity for `_atom_site` loops so that noisy metadata loops do not abort a read.
- **Typed atom storage**: `atom::AtomSite` accumulates `_atom_site.*` values with numeric parsing that understands estimated standard deviations (`1.234(5)`). The helper returns `AtomPosition` for downstream geometry, similar to `pdbtbx::Atom` but without mutation-heavy APIs.
- **Unit cell builder**: `_cell` values are buffered until all six parameters are available. Partial definitions are rejected with contextful `ParseError`, while fully missing cells simply yield `None`.
- **Inter-atom analysis**: `model::MmcifStructure::inspect_contacts` walks pairwise atoms, optionally applying minimum-image wrapping for orthogonal cells (`UnitCell::minimum_image_distance`). Contacts are classified using the in-crate `ElementData` table, giving quick covalent/clash diagnostics without pulling in a force-field library.
- **Async + gzip support**: `parse_file`, `parse_reader`, `parse_async_reader` share the same parser core. Files ending with `.gz` automatically stream through `flate2` so `.cif.gz` bundles (like `11AS.cif.gz` shipped at repo root) can be read transparently.

Tests cover both synthetic strings and the real `11AS.cif.gz` file to ensure the parser handles field noise, compression, and realistic atom counts.
