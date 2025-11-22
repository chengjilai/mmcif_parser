import pathlib
import unittest

import mmcif_parser

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[3]
CIF_PATH = PROJECT_ROOT / "11AS.cif.gz"


class StructureUtilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if not CIF_PATH.is_file():
            raise unittest.SkipTest(f"test data missing: {CIF_PATH}")
        cls.structure = mmcif_parser.parse(CIF_PATH)

    def test_unit_cell_metadata(self) -> None:
        unit_cell = self.structure.unit_cell()
        self.assertIsNotNone(unit_cell)
        lengths, angles = unit_cell
        self.assertEqual(len(lengths), 3)
        self.assertEqual(len(angles), 3)
        self.assertTrue(all(value > 0.0 for value in lengths))
        self.assertTrue(all(0.0 < value < 180.0 for value in angles))

    def test_contacts_sorted_and_classified(self) -> None:
        contacts = self.structure.contacts(2.2, max_results=5)
        self.assertTrue(contacts)
        self.assertLessEqual(len(contacts), 5)
        distances = [contact.distance for contact in contacts]
        self.assertEqual(distances, sorted(distances))
        valid_classes = {"covalent", "short_contact", "clash", "nonbonded"}
        self.assertTrue(all(contact.classification in valid_classes for contact in contacts))

    def test_atom_profile_neighbors_and_properties(self) -> None:
        profile = self.structure.atom_profile(0, include_geometric=True, tolerance=0.4)
        self.assertEqual(profile.index, 0)
        self.assertIsNotNone(profile.label_comp_id)
        self.assertTrue(profile.neighbor_indices)
        self.assertEqual(len(profile.neighbor_indices), len(profile.neighbor_bonds))
        self.assertIn(
            profile.hybridization,
            {None, "sp", "sp2", "sp3", "aromatic"},
        )


if __name__ == "__main__":
    unittest.main()
