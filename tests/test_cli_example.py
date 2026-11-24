"""End-to-end run of the bundled MIC example."""
import tempfile
import unittest
from pathlib import Path

from _support import EXAMPLES, read_csv, run_cli, run_folder


class ExampleRunTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._temp = tempfile.TemporaryDirectory()
        work = Path(cls._temp.name)
        code, out, err = run_cli("--workbook", EXAMPLES / "demo_mic.xlsx", "--output", work / "out", cwd=work)
        if code:
            raise AssertionError(f"CLI exited {code}\n{err[-2000:]}")
        cls.run_dir = run_folder(out)

    @classmethod
    def tearDownClass(cls):
        cls._temp.cleanup()

    def test_agar_dilution_mics(self):
        from openpyxl import load_workbook

        table = load_workbook(self.run_dir / "01_Mutation_Analysis/NG_AMR_MIC.xlsx", read_only=True)
        try:
            header, *rows = table.active.values
            mic = header.index("Penicillin G | AD MIC (µg/mL)")
            self.assertEqual({r[0]: r[mic] for r in rows}, {"DEMO-A": 0.25, "DEMO-B": 0.5, "DEMO-C": 1})
        finally:
            table.close()

    def test_expected_outputs_exist(self):
        for rel in ("01_Mutation_Analysis/NG_AMR_MIC.xlsx",
                    "05_Other_Outputs/NG_AMR_Run_Manifest.json",
                    "05_Other_Outputs/run_status.csv",
                    "05_Other_Outputs/analysis.log"):
            self.assertTrue((self.run_dir / rel).is_file(), rel)

    def test_no_images_or_slides(self):
        made = [p.name for p in self.run_dir.rglob("*") if p.suffix.lower() in {".png", ".pptx", ".docx"}]
        self.assertEqual(made, [])

    def test_workbook_sheets(self):
        from openpyxl import load_workbook

        book = load_workbook(self.run_dir / "01_Mutation_Analysis/NG_AMR_MIC.xlsx", read_only=True)
        try:
            self.assertEqual(book.sheetnames, ["Isolate MICs", "Method details", "ID reconciliation", "QC",
                                               "Source cells", "Paired agreement", "Recorded WGS agreement",
                                               "Method pairs", "AD comparator errors"])
        finally:
            book.close()

    def test_no_failed_stages(self):
        rows = read_csv(self.run_dir / "05_Other_Outputs/run_status.csv")
        self.assertTrue(rows)
        self.assertFalse([r for r in rows if r["Status"] == "failed"])


if __name__ == "__main__":
    unittest.main()
