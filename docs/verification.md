# Verification

The source release passed the following Linux/Python 3.12 checks with fictional fixtures.

| Workflow | Check | Result |
|---|---|---|
| AMR | Read the five-sheet workbook, verify three named isolates and their AD MIC values, and generate the comparison and summary presentations. | Passed |
| AMR regression checks | Run the inherited self-test with Biopython available. | All 33 checks passed |

To run the example test in this repository, install its requirements, then run `python -m unittest discover -s tests -v`. These checks cover software behavior. The source release did not validate clinical performance, full assembly calling with BLAST+, or Docker-based tree inference.
