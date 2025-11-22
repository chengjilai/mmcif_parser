#!/usr/bin/env python3
# type: ignore
# pyright: reportGeneralTypeIssues=false
"""Summarize an mmCIF file using Biopython for cross-checks."""

from __future__ import annotations

import gzip
import json
import pathlib
import sys
import types
from typing import Any, cast

if len(sys.argv) != 2:
    raise SystemExit("usage: biopython_summary.py <path-to-cif.gz>")

repo_root = pathlib.Path(__file__).resolve().parents[3]
source_biopython = repo_root / "biopython"


def _ensure_cealign_stub() -> None:
    ccealign_module = "Bio.PDB.ccealign"
    if ccealign_module in sys.modules:
        return
    stub = types.ModuleType(ccealign_module)

    def _missing(*_args: Any, **_kwargs: Any) -> None:  # pragma: no cover - shim
        raise ImportError("Bio.PDB.ccealign extension is not available")

    stub.run_cealign = _missing  # type: ignore[attr-defined]
    sys.modules[ccealign_module] = stub


def _load_biopython():
    try:
        from Bio.PDB.MMCIFParser import MMCIFParser  # type: ignore
        from Bio.PDB.Model import Model  # type: ignore
        from Bio.PDB.Structure import Structure  # type: ignore
        return MMCIFParser, Model, Structure
    except ModuleNotFoundError:
        if not source_biopython.exists():
            raise
        sys.path.insert(0, str(source_biopython))
        _ensure_cealign_stub()
        from Bio.PDB.MMCIFParser import MMCIFParser  # type: ignore
        from Bio.PDB.Model import Model  # type: ignore
        from Bio.PDB.Structure import Structure  # type: ignore
        return MMCIFParser, Model, Structure


try:
    MMCIFParser, Model, Structure = _load_biopython()
except ModuleNotFoundError as exc:  # pragma: no cover - clear error for Rust test
    raise SystemExit(
        "Biopython is required for this test (pip install biopython)."
    ) from exc

parser = MMCIFParser(QUIET=True)
with gzip.open(sys.argv[1], "rt") as handle:
    structure = cast(Structure, parser.get_structure("reference", handle))

model = cast(Model, next(structure.get_models()))
atoms = list(model.get_atoms())
residues = list(model.get_residues())
chains = sorted({chain.id.strip() for chain in model.get_chains() if chain.id.strip()})
residue_names = sorted(
    {
        residue.get_resname().strip().upper()
        for residue in residues
        if residue.get_resname().strip()
    }
)

payload = {
    "atom_count": len(atoms),
    "residue_count": len(residues),
    "chain_ids": chains,
    "residue_names_sample": residue_names[:10],
    "first_atom_name": atoms[0].get_name().strip(),
    "first_residue": residues[0].get_resname().strip().upper(),
}

print(json.dumps(payload))
