# mmcif-parser (Python bindings)

Python wheels generated from the `mmcif-core` and `mmcif-analyze` Rust crates. The
extension module exposes a high-level `Structure` type that mirrors the Rust data
model and offers convenience helpers for spatial queries.

```python
from pathlib import Path
import mmcif-parser as mmcif

structure = mmcif.parse(Path("11AS.cif.gz"))
print(structure.atom_count)
for hit in structure.nearest_atoms((12.0, 4.2, 0.3), max_distance=5.0, max_results=5):
    print(hit.index, hit.distance)
```

## Building

1. Create a virtual environment and install maturin once:

    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    pip install maturin
    ```

2. Use the provided automation targets:

    ```bash
    cd python/mmcif-parser_py
    make develop   # editable install via maturin develop --release
    make build     # produce wheels + sdist via maturin build --release --sdist --strip
    ```

`make distclean` removes the build artifacts and uninstalls the editable package
from the virtual environment, while `make clean` only clears the local `target`
and `egg-info` directories.
