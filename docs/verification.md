# Verification

Run the checks with:

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
.\.venv\Scripts\python.exe run_amr.py --self-test
```

| Check | What it covers | Data |
|---|---|---|
| `test_cli_example` | Reads the five-sheet example workbook, checks agar-dilution MICs of three isolates, workbook sheets, that no stage failed and that no image or slide files are made | Fictional workbook |
| `test_cli_args` | Rejects missing inputs, missing files and invalid thresholds | — |
| `test_units` | QC-strain name variants, alias file rules, MIC value parsing, workbook reading | Fictional |
| `test_selftest` | The 23 built-in checks: protein and nucleotide mutation calling, disruption rules, gene types, prevalence threshold arithmetic, logger formats | Synthetic sequences |
| `test_genomic_demo` | Full BLAST+ route: unchanged isolates have no calls; penA A10V and gyrA S91F are called and flagged as curated; a premature stop is reported as disruption; a missing locus is reported as absent | Random synthetic genomes with planted changes (`tests/fixtures/make_genomic_demo.py`) |

When this repository was restructured (v1.0.0), the tables, CSVs, PDFs and HTML reports of a full example
run were compared file by file against the original code's output, ignoring timestamps.

## Not validated

These checks cover software behaviour on fictional and synthetic inputs. They do **not** establish clinical
or diagnostic performance, accuracy on real *N. gonorrhoeae* assemblies, or the correctness of any curated
variant catalogue you supply. Review results before using them.
