"""Run the built-in regression checks under unittest."""
import logging
import unittest

import _support  # noqa: F401  (puts src/ on sys.path)

from ng_amr.genomics.alignment import Seq
from ng_amr.selftest import run_self_tests


class SelfTest(unittest.TestCase):
    def test_built_in_checks(self):
        logging.disable(logging.CRITICAL)
        try:
            verdict = run_self_tests()
        finally:
            logging.disable(logging.NOTSET)
        if Seq is None:
            self.assertIsNone(verdict, "core checks failed")
            self.skipTest("Biopython not installed: alignment-dependent checks skipped")
        self.assertTrue(verdict)


if __name__ == "__main__":
    unittest.main()
