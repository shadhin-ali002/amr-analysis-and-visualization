"""Small unit tests of identifier handling and workbook parsing."""
import tempfile
import unittest
from pathlib import Path

from _support import EXAMPLES

from ng_amr.core import ids
from ng_amr.core.audit import Audit
from ng_amr.core.constants import QC
from ng_amr.io.workbook import measurement, read_master


class IdentifierTest(unittest.TestCase):
    def tearDown(self):
        ids.ALIASES.clear()

    def test_qc_strain_spellings(self):
        for raw in ("ATCC49226", "atcc 49226", "QC-ATCC_49226", "ATCC49226 control"):
            with self.subTest(raw=raw):
                self.assertEqual(ids.canonical(raw), QC)

    def test_canonical_upper_cases_and_strips(self):
        self.assertEqual(ids.canonical("  demo-a "), "DEMO-A")

    def test_alias_file(self):
        aliases = ids.load_aliases(EXAMPLES / "demo_aliases.csv")
        self.assertEqual(aliases, {"DEMO-OLD": "DEMO-A"})
        ids.ALIASES.update(aliases)
        self.assertEqual(ids.canonical("demo-old"), "DEMO-A")

    def test_alias_file_rejects_chains_and_qc(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "aliases.csv"
            for body in ("alias,canonical_id\nA,B\nB,C\n", "alias,canonical_id\nX,ATCC49226\n",
                         "wrong,header\nA,B\n"):
                with self.subTest(body=body):
                    path.write_text(body, encoding="utf-8")
                    with self.assertRaises(ValueError):
                        ids.load_aliases(path)

    def test_assembly_id_drops_numeric_prefix(self):
        self.assertEqual(ids.assembly_id("12_demo-b"), "DEMO-B")


class WorkbookTest(unittest.TestCase):
    def test_measurement(self):
        self.assertEqual(measurement("0.25"), ("", 0.25))
        self.assertEqual(measurement("<=0.016"), ("<=", 0.016))
        self.assertEqual(measurement(">256"), (">", 256.0))
        self.assertEqual(measurement("1e-1"), ("", 0.1))
        for bad in ("", "n/a", "0", "-1"):
            with self.subTest(bad=bad):
                self.assertEqual(measurement(bad), ("", None))

    def test_read_demo_workbook(self):
        audit = Audit()
        book = read_master(EXAMPLES / "demo_mic.xlsx", audit)
        self.assertFalse([row for row in audit.rows if row[1] == "failed"])
        self.assertTrue({"DEMO-A", "DEMO-B", "DEMO-C"} <= set().union(*book.ids.values()))


if __name__ == "__main__":
    unittest.main()
