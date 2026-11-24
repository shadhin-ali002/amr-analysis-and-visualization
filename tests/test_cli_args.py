"""Argument validation of the command line."""
import tempfile
import unittest
from pathlib import Path

from _support import EXAMPLES, run_cli


class ArgumentTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.work = Path(temp.name)

    def assert_usage_error(self, *args, message):
        code, _, err = run_cli(*args, "--output", self.work / "out", cwd=self.work)
        self.assertEqual(code, 2)
        self.assertIn(message, err)
        self.assertFalse((self.work / "out").exists())

    def test_requires_an_input(self):
        self.assert_usage_error(message="Supply --workbook or --wgs-base")

    def test_missing_workbook(self):
        self.assert_usage_error("--workbook", self.work / "absent.xlsx", message="Missing input file")

    def test_bad_thresholds(self):
        for value in ("0,50", "abc", "101"):
            with self.subTest(value=value):
                self.assert_usage_error("--workbook", EXAMPLES / "demo_mic.xlsx", "--thresholds", value,
                                        message="Thresholds must be whole percentages")

    def test_version(self):
        from ng_amr import __version__

        code, out, _ = run_cli("--version", cwd=self.work)
        self.assertEqual(code, 0)
        self.assertIn(__version__, out)


if __name__ == "__main__":
    unittest.main()
