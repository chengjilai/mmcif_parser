# pyright: reportMissingImports=false, reportMissingModuleSource=false, reportGeneralTypeIssues=false

import gzip
import math
import pathlib
import sys
import unittest
from importlib import import_module
from typing import Any

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[3]
CIF_PATH = PROJECT_ROOT / "11AS.cif.gz"
BIOPYTHON_PATH = PROJECT_ROOT / "biopython"
if BIOPYTHON_PATH.exists():
    sys.path.insert(0, str(BIOPYTHON_PATH))

try:
    _bio_pdb = import_module("Bio.PDB")
except ModuleNotFoundError:  # pragma: no cover - best-effort import
    BioMMCIFParserCls: Any = None
else:
    BioMMCIFParserCls = getattr(_bio_pdb, "MMCIFParser", None)


def _atom_coord(atom: Any) -> tuple[float, float, float]:
    coord = atom.coord  # type: ignore[attr-defined]
    return (float(coord[0]), float(coord[1]), float(coord[2]))


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(
        (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
    )


@unittest.skipIf(BioMMCIFParserCls is None, "Biopython MMCIFParser unavailable")
class BiopythonParityTests(unittest.TestCase):
    structure: Any
    bio_structure: Any
    first_residue: Any
    ca_atom: Any
    n_atom: Any

    @classmethod
    def setUpClass(cls) -> None:
        if not CIF_PATH.is_file():
            raise unittest.SkipTest(f"test data missing: {CIF_PATH}")
        mmcif_module = import_module("mmcif-parser")
        cls.structure = mmcif_module.parse(CIF_PATH)  # type: ignore[attr-defined]
        if BioMMCIFParserCls is None:
            raise unittest.SkipTest("Biopython MMCIFParser unavailable")
        parser = BioMMCIFParserCls(QUIET=True)
        with gzip.open(CIF_PATH, "rt") as handle:
            cls.bio_structure = parser.get_structure("11AS", handle)
        cls.first_residue = next(cls.bio_structure.get_residues())
        cls.ca_atom = cls.first_residue["CA"]
        cls.n_atom = cls.first_residue["N"]

    def test_atom_count_matches_biopython(self) -> None:
        bio_atom_count = sum(1 for _ in self.bio_structure.get_atoms())
        self.assertEqual(self.structure.atom_count, bio_atom_count)

    def test_nearest_atom_matches_biopython_ca(self) -> None:
        ca_origin = _atom_coord(self.ca_atom)
        hits = self.structure.nearest_atoms(ca_origin, max_distance=0.5, max_results=1)
        self.assertEqual(len(hits), 1)
        hit = hits[0]
        self.assertEqual(hit.label_atom_id, self.ca_atom.get_name().strip())
        self.assertEqual(hit.label_comp_id, self.first_residue.get_resname().strip())
        self.assertEqual(hit.chain_id, self.ca_atom.get_parent().get_parent().id)
        self.assertEqual(hit.seq_id, int(self.first_residue.get_id()[1]))
        self.assertAlmostEqual(hit.distance, 0.0, places=5)

    def test_distance_matrix_matches_pair_distance(self) -> None:
        ca_origin = _atom_coord(self.ca_atom)
        n_origin = _atom_coord(self.n_atom)
        ca_hit = self.structure.nearest_atoms(ca_origin, max_distance=0.5, max_results=1)[0]
        n_hit = self.structure.nearest_atoms(n_origin, max_distance=0.5, max_results=1)[0]
        matrix = self.structure.dense_distance_matrix([ca_hit.index, n_hit.index])
        self.assertEqual(matrix.atom_indices, [ca_hit.index, n_hit.index])
        observed = matrix.distances[0][1]
        expected = _distance(ca_origin, n_origin)
        self.assertAlmostEqual(observed, expected, places=5)
        bonds = self.structure.explicit_bonds()
        normalized = {(min(i, j), max(i, j)) for i, j in bonds}
        key = (min(ca_hit.index, n_hit.index), max(ca_hit.index, n_hit.index))
        self.assertIn(key, normalized)


if __name__ == "__main__":
    unittest.main()
