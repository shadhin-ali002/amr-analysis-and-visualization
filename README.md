# AMR analysis and visualization

[![Tests](https://github.com/shadhin-ali002/amr-analysis-and-visualization/actions/workflows/tests.yml/badge.svg)](https://github.com/shadhin-ali002/amr-analysis-and-visualization/actions/workflows/tests.yml)
![Python](https://img.shields.io/badge/python-3.11%2B-blue)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

`ng-amr` is a command-line tool for antimicrobial-resistance (AMR) analysis of *Neisseria gonorrhoeae*
isolates. It does two things:

1. **MIC workbook analysis.** It reads a laboratory MIC workbook (E-strip, agar dilution, McFarland,
   disc diffusion and recorded WGS categories) and writes one clean Excel workbook. That workbook has
   per-isolate MICs, QC results, source-cell traceability, E-strip vs agar-dilution agreement (essential and
   categorical agreement, log₂ bias, Wilcoxon test) and method-pair error rates.
2. **Resistance-gene mutation detection** (optional). From existing genome assemblies it locates 30
   resistance loci with BLAST+ (penA, mtrR, gyrA, parC, 23S rDNA and others), aligns them to reference
   sequences, and calls substitutions, insertions, deletions, premature stops and promoter or rRNA SNPs. Calls
   are flagged against a curated variant catalogue, then written as Excel matrices, CSV tables, and per-gene
   and per-isolate **PDF and HTML alignment reports**.

The tool makes no image files or slide decks. Every result is a table (XLSX/CSV), a report (PDF/HTML) or
plain text (FASTA/Markdown).

## Install

Python 3.11 or newer. On Windows PowerShell:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[genomics]"
```

Use `pip install -e .` if you only need the MIC workbook analysis. The genomic route also needs
[NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (`blastn`,
`makeblastdb`) on `PATH`. On macOS or Linux use `python3 -m venv .venv` and `.venv/bin/python`.

Without installing, you can also run `pip install -r requirements-genomics.txt` and call `python run_amr.py`
instead of `ng-amr`.

## Quick start

```powershell
ng-amr --workbook examples/demo_mic.xlsx
ng-amr --workbook examples/demo_mic.xlsx --wgs-base tests/fixtures/genomic_demo
ng-amr --self-test
```

The first command analyses the fictional example workbook. The second also runs mutation detection on a
small synthetic genome set, which takes about 15 seconds. `ng-amr --help` lists every option.

| Option | Purpose |
|---|---|
| `--workbook` | MIC master workbook (`.xlsx`) |
| `--wgs-base` | Folder holding `Assembly/NT/Plane_Fasta/*.fasta` and `Refs/` |
| `--catalog` | Curated variant JSON (default `Refs/Curated_Variants.json`) |
| `--genes` | Comma-separated subset of loci, e.g. `penA,gyrA` |
| `--id-map` | `filename,isolate_id` CSV when FASTA names differ from isolate IDs |
| `--aliases` | `alias,canonical_id` CSV to reconcile historical IDs |
| `--thresholds` | Prevalence percentages for summary tables (default `25,30,35,40,50`) |
| `--no-html` | Skip HTML alignment reports |
| `--output` | Parent folder for timestamped runs (default `outputs/`) |

## Outputs

Each run writes a new `outputs/Run_<timestamp>/` folder:

```
Run_<timestamp>/
├── 01_Mutation_Analysis/
│   ├── NG_AMR_MIC.xlsx                     MIC tables, QC, agreement statistics
│   ├── NG_AMR_Mutation_MIC_Matrix.xlsx     mutation × isolate matrix with MICs   (genomic)
│   ├── NG_AMR_MIC_MasterMatrix.xlsx        combined master workbook              (genomic)
│   ├── WGS_calls.csv                       locus status per isolate              (genomic)
│   ├── Mutation_Tables/                    per-gene call and matrix CSVs         (genomic)
│   ├── Mutation_Reports/By_Gene|By_Isolate PDF, HTML and Excel alignment reports  (genomic)
│   ├── Extracted_Sequences/                called gene and protein FASTA         (genomic)
│   └── Summary_Tables/                     catalogue, burden, thresholds, drug classes, key findings
└── 05_Other_Outputs/
    ├── run_status.csv                      outcome of every stage
    ├── NG_AMR_Run_Manifest.json            arguments, versions, input and source checksums
    └── analysis.log
```

See [docs/outputs.md](docs/outputs.md) for every file and [docs/data-formats.md](docs/data-formats.md) for the inputs.

## Project layout

```
src/ng_amr/
├── cli.py               command line and run sequence
├── core/                identifiers, aliases, audit trail, gene and antibiotic configuration
├── io/                  MIC workbook reader, assembly discovery
├── phenotype/           MIC workbook and agreement statistics
├── genomics/            BLAST+ search, alignment, variant calling, sequence export
├── reports/             PDF / HTML / Excel mutation reports, master matrices, summary tables
└── selftest.py          built-in regression checks
tests/                   unit, CLI and end-to-end tests (incl. synthetic BLAST+ fixture)
examples/                fictional demo workbook and alias file
```

## Testing

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

The suite runs the example end to end and checks the CLI's argument validation, identifier and workbook
parsing, and the 23 built-in checks. It also runs mutation detection on the synthetic genome set, which
contains planted changes (penA A10V, gyrA S91F, a premature stop, a missing locus, promoter and 23S SNPs).
The BLAST+ test is skipped when BLAST+ or the genomics extras are not installed.
[docs/verification.md](docs/verification.md) says what is and is not validated.

## Data privacy

All example and test data are fictional or randomly generated. Keep real isolate data in `data/` or
`private/`. Both folders, and `outputs/`, are git-ignored.

## Citation

If you use this software, please cite it using [CITATION.cff](CITATION.cff).

## License

[MIT](LICENSE) © 2026 Md. Shadhin Ali
