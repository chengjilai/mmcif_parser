# mmcif_parser

Rust crates and Python bindings for fast mmCIF parsing, structural analysis, and
high-level workflows. This repo now keeps a single top-level README; each module’s
original documentation remains in place and is linked below for deeper dives.

## Repository Layout

| Path | Purpose | Linked README |
| --- | --- | --- |
| `crates/mmcif_core` | Core mmCIF reader, data model, and utilities | [Module README](crates/mmcif_core/README.md) |
| `crates/mmcif_analyze` | Spatial search, bond inference, and chemical profiling on top of `mmcif_core` | [Module README](crates/mmcif_analyze/README.md) |
| `python/mmcif_parser_py` | PyO3 bindings packaged as the `mmcif_parser` wheel | [Module README](python/mmcif_parser_py/README.md) |
| `docs/` | Documentation guide plus generated artifacts | [Docs README](docs/README.md) |

All other README files should link back to this page where appropriate to keep
navigation consistent.

## Quick Start

1. **Rust toolchain**: Install Rust ≥ 1.85 via [rustup](https://rustup.rs/), then
   build/test every crate with:
   ```bash
   cargo fmt
   cargo test --workspace
   ```
2. **Python bindings**: Create a virtual environment and install via `maturin`:
   ```bash
   cd python/mmcif_parser_py
   python3 -m venv .venv && source .venv/bin/activate
   pip install maturin
   make develop
   python -c "import mmcif_parser; print(mmca := mmcif_parser.parse('11AS.cif.gz').atom_count)"
   ```
3. **Documentation**: Generate Rust + optional Python API docs using the helper
   script described in [docs/README.md](docs/README.md):
   ```bash
   make docs
   ```

## Release Checklist

- Update crate versions and the Python wheel metadata as needed.
- Run `cargo test --workspace` and `pytest` (or `python -m unittest`) inside
  `python/mmcif_parser_py/tests`.
- Regenerate docs (`make docs`) if public APIs changed.
- Publish Rust crates (if desired) and the Python wheel (`maturin publish`).

## Contributing

Issues and pull requests are welcome. Please ensure new functionality includes
unit tests (Rust or Python as appropriate) and update the relevant module README
so that this top-level overview remains the authoritative entry point.
