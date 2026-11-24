"""Genomic route on the synthetic fixture: BLAST+, alignment and variant calls.

Skipped unless Biopython, parasail and NCBI BLAST+ are installed. The fixture is random
sequence with planted changes (see fixtures/make_genomic_demo.py).
"""
import importlib.util
import shutil
import tempfile
import unittest
from pathlib import Path

from _support import EXAMPLES, FIXTURES, read_csv, run_cli, run_folder

READY = (all(importlib.util.find_spec(m) for m in ("Bio", "parasail"))
         and shutil.which("blastn") and shutil.which("makeblastdb"))


@unittest.skipUnless(READY, "needs biopython, parasail and BLAST+ on PATH")
class GenomicDemoTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._temp = tempfile.TemporaryDirectory()
        work = Path(cls._temp.name)
        code, out, err = run_cli("--workbook", EXAMPLES / "demo_mic.xlsx",
                                 "--wgs-base", FIXTURES / "genomic_demo",
                                 "--genes", "penA,gyrA", "--no-html", "--output", work / "out", cwd=work)
        if code:
            raise AssertionError(f"CLI exited {code}\n{err[-3000:]}")
        cls.run_dir = run_folder(out)
        tables = cls.run_dir / "01_Mutation_Analysis" / "Mutation_Tables"
        cls.calls = {gene: {r["Isolate"]: r for r in read_csv(tables / f"{gene}_protein_summary.csv")}
                     for gene in ("penA", "gyrA")}

    @classmethod
    def tearDownClass(cls):
        cls._temp.cleanup()

    def test_reference_identical_isolates_have_no_mutations(self):
        for gene in ("penA", "gyrA"):
            for isolate in ("DEMO-A", "ATCC49226"):
                with self.subTest(gene=gene, isolate=isolate):
                    row = self.calls[gene][isolate]
                    self.assertEqual(row["Status"], "intact")
                    self.assertEqual(row["Detail"], "No mutations")

    def test_curated_substitutions(self):
        for gene, mutation in (("penA", "A10V"), ("gyrA", "S91F")):
            with self.subTest(gene=gene):
                row = self.calls[gene]["DEMO-B"]
                self.assertEqual((row["Mutation"], row["Type"], row["Curated"]),
                                 (mutation, "substitution", "Curated"))

    def test_premature_stop_is_disruption(self):
        row = self.calls["penA"]["DEMO-C"]
        self.assertEqual(row["Status"], "disrupted")
        self.assertIn("premature stop at aa 120", row["Detail"])

    def test_absent_locus(self):
        self.assertEqual(self.calls["gyrA"]["DEMO-C"]["Status"], "absent")

    def test_no_failed_stages(self):
        rows = read_csv(self.run_dir / "05_Other_Outputs/run_status.csv")
        self.assertFalse([r for r in rows if r["Status"] == "failed"])


if __name__ == "__main__":
    unittest.main()
