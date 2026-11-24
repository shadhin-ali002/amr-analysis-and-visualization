"""Exercise the project CLI on fictional input in a separate working directory."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


class ExampleTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.work = Path(temp.name)

    def run_workflow(self, *args):
        result = subprocess.run([sys.executable, ROOT / "run_amr.py", *args, "--output", self.work / "out"],
                                cwd=self.work, capture_output=True, text=True)
        if result.returncode:
            raise subprocess.CalledProcessError(result.returncode, result.args,
                                                output=result.stdout, stderr=result.stderr)
        return Path(result.stdout.strip().splitlines()[-1])

    def test_example(self):
        from openpyxl import load_workbook
        run = self.run_workflow("--workbook", EX / "demo_mic.xlsx")
        table = load_workbook(run / "01_Mutation_Analysis/NG_AMR_MIC.xlsx", read_only=True)
        try:
            header, *rows = table.active.values
            mic = header.index("Penicillin G | AD MIC (µg/mL)")
            self.assertEqual({r[0]:r[mic] for r in rows}, {"DEMO-A":0.25, "DEMO-B":0.5, "DEMO-C":1})
        finally:
            table.close()
        self.assertTrue(list(run.rglob("*.pptx")))


if __name__ == "__main__":
    unittest.main()
