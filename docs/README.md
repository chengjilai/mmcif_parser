# Documentation Guide

> For the full project overview, see the [root README](../README.md). This page
> only covers documentation tooling.

This repository ships both Rust and Python components. To keep the generated
API references consistent, run the helper script under `scripts/`:

```bash
./scripts/generate_docs.sh
```

The script performs the following steps:

1. Builds the Rust HTML documentation for every crate in the workspace via
   `cargo doc --workspace --no-deps`. The output lives under `target/doc`.
2. If the `mmcif_parser` Python module is importable and `pdoc` is installed in
   your current environment, it emits Markdown/HTML documentation into
   `docs/api/python`. When either prerequisite is missing, the script logs a
   warning and skips the Python step so that `cargo doc` is still available.

## Prerequisites

- A recent stable Rust toolchain (Rust 1.80+). Install with
  [`rustup`](https://rustup.rs/) if you have not already.
- (Optional) `pdoc` ≥ 14.0 in the Python environment from which you run the
  script. Install it via `pip install pdoc`. The module `mmcif_parser` must also
  be importable (for example by running `make develop` inside
  `python/mmcif_parser_py`).

## Cleaning generated files

Rust docs are stored under `target/doc` (already git-ignored). Python docs are
placed under `docs/api/python`. Remove the latter with `rm -rf docs/api/python`
if you need a clean slate or run `make clean-docs` from the repository root once
that target is available in your environment.
