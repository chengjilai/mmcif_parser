#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY_OUTPUT_DIR="$ROOT_DIR/docs/api/python"
PYTHON_BIN="${PYTHON_BIN:-python3}"

pushd "$ROOT_DIR" >/dev/null

echo "[docs] Building Rust API reference via cargo doc" >&2
cargo doc --workspace --no-deps

mkdir -p "$PY_OUTPUT_DIR"

if "$PYTHON_BIN" -m pdoc --version >/dev/null 2>&1; then
    if "$PYTHON_BIN" - <<'PY' >/dev/null 2>&1; then
import importlib.util
import sys
spec = importlib.util.find_spec("mmcif_parser")
sys.exit(0 if spec is not None else 1)
PY
        then
            echo "[docs] Generating Python API documentation via pdoc" >&2
            "$PYTHON_BIN" -m pdoc mmcif_parser --output-dir "$PY_OUTPUT_DIR"
        else
            echo "[docs] mmcif_parser is not importable; skipping Python docs" >&2
        fi
    else
        echo "[docs] Failed to probe pdoc availability; skipping Python docs" >&2
    fi
else
    echo "[docs] pdoc is not installed (try 'pip install pdoc'); skipping Python docs" >&2
fi

popd >/dev/null
