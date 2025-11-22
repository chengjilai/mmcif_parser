Changelog of mmcif_parser
document every small performance progress

## Unreleased

- `BondEdge` is now a `Copy` type reused directly by sparse adjacencies, so both the `mmcif_core` cache and `mmcif_analyze` avoid intermediate cloning when wiring neighbors.
- `compute_bonds_by_distance` precomputes the position/radius tuple for each usable atom and hoists the unit-cell/Euclidean branch outside the nested loop, reducing per-pair work in dense structures.
- The sparse adjacency builder starts from the typed explicit cache and only appends geometric edges when they are missing, preventing duplicate entries while keeping the fast path for distance-derived bonds.
- Criterion bench harnesses for each crate are registered via `[[bench]]` sections, allowing direct `cargo bench --bench …` execution without the default test harness.
- `nearest_amino_acids` preallocates its residue map and matches on the unit cell once, cutting per-atom option checks while keeping the minimal-distance update logic.
- `build_adjacency` now performs a lightweight degree count to reserve the exact neighbor capacity for each atom before inserting edges, eliminating repeated reallocations during dense builds.
- The Criterion suite tracks dense adjacency construction and nearest-residue queries in addition to the existing sparse, dense-distance, and nearest-atom runs.
- `nearest_amino_acids` now caches each atom’s residue key before distance evaluation, so repeated inserts reuse the derived identifiers instead of re-reading the labeling tags.
- Added `python/mmcif_parser_py`, a `pyo3` + `maturin` binding crate that exposes parsing, spatial queries, and adjacency helpers as a PyPI-ready `mmcif_parser` module.
- Added a `Makefile` in `python/mmcif_parser_py` with `make develop`, `make build`, `make clean`, and `make distclean` targets to automate maturin workflows and venv cleanup.
- Python bindings now export unit-cell metadata, contact inspection, and per-atom property profiles with accompanying regression tests built on the `11AS.cif.gz` structure.
- Added Biopython-backed parity tests (`test_biopython_parity.py`) to cross-check atom counts, residue metadata, and pairwise distances against the reference MMCIF parser using the bundled 11AS structure.