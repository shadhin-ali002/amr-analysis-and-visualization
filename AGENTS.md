# AGENTS.md

Guidance for AI coding agents working in this repository.

## What this is

`amr-analysis-and-visualization` (package `ng_amr`, command `ng-amr`) is a Python command-line tool for
antimicrobial-resistance (AMR) analysis of *Neisseria gonorrhoeae* isolates, published by Md. Shadhin Ali
(`shadhin-ali002`). It:

- reads an MIC workbook (sheets `EST_MIC`, `AD_MIC`, `McFarland`, `Disc_Diffusion`, `WGS`) and writes
  `NG_AMR_MIC.xlsx` with MIC tables, QC, and method-agreement statistics;
- optionally, with genome assemblies and reference files, finds 30 resistance loci with BLAST+, calls
  variants, and writes Excel matrices, CSV tables and PDF/HTML alignment reports.

**Outputs are tables (XLSX/CSV), reports (PDF/HTML) and text only.** Image, PowerPoint and Word generation
was removed on purpose in v1.0.0. Do not add matplotlib/PNG/PPTX/DOCX output back unless the owner asks.

It is one of five sibling repos (`e-strip-plate-analysis`, `agar-dilution-plate-analysis`,
`phylogeny-analysis`, `virulence-analysis`). Each must run on its own. The owner's profile README and
website live in `../shadhin-ali002` and `../shadhin-ali002.github.io`.

## Files

| Path | Purpose |
| --- | --- |
| `src/ng_amr/cli.py` | Argument parsing and the run sequence (`ng-amr`, `python -m ng_amr`) |
| `src/ng_amr/core/` | `ids.py` (canonical IDs, aliases, QC strain ATCC 49226), `audit.py` (run audit, logger, CSV helper, `stage`), `constants.py` (antibiotics, gene panel, gene types) |
| `src/ng_amr/io/` | `workbook.py` (MIC workbook reader), `assemblies.py` (FASTA discovery) |
| `src/ng_amr/phenotype/mic.py` | `NG_AMR_MIC.xlsx`: MIC sheets, paired/categorical agreement, method-pair errors |
| `src/ng_amr/genomics/` | `locus.py` (BLAST+, per-locus calling, `process_gene`), `alignment.py` (translation, alignment, mutation calling), `sequences.py`, `aa_comparison.py`, `pipeline.py` (genomic stage order) |
| `src/ng_amr/reports/` | `pdf.py`, `html.py`, `grouped.py` (by-gene / by-isolate reports), `master_matrix.py`, `summary_tables.py`, `run_summary.py` (manifest) |
| `src/ng_amr/selftest.py` | 23 built-in checks (`ng-amr --self-test`) |
| `run_amr.py` | Runs the CLI from a checkout without installing |
| `tests/` | `unittest` suite; `fixtures/genomic_demo` is a random synthetic genome set made by `fixtures/make_genomic_demo.py` |
| `examples/` | Fictional demo workbook (`DEMO-A/B/C`, penicillin G) and alias CSV |
| `docs/` | `data-formats.md`, `outputs.md`, `verification.md` |
| `pyproject.toml` | Metadata, dependencies (`[genomics]` extra), ruff config |
| `LICENSE` | MIT. Don't change it, or add code under an incompatible licence, without the owner's request |

## Setup and checks

Windows PowerShell, Python 3.11+:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[genomics]" ruff
.\.venv\Scripts\python.exe -m ruff check src tests
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
.\.venv\Scripts\python.exe run_amr.py --self-test
```

Expected: ruff clean, all tests pass (about 10 s; the BLAST+ test is skipped when BLAST+ is not on PATH),
and `ALL SELF-TESTS PASSED (23 checks)`. CI (`.github/workflows/tests.yml`) runs the same on Ubuntu and
Windows. Write any manual run output outside the repo or to `outputs/` (git-ignored).

## Rules

- **Do not change results.** This produces scientific numbers. Refactors must keep every output file name,
  folder name, column header and value identical. Compare the example output before and after a change.
  Changes to calling logic, thresholds, breakpoints or the curated catalogue need the owner's explicit request.
- **Move code, don't rewrite it.** When restructuring, cut and paste functions unchanged, then re-run the checks.
  Keep `cli.py` imports lazy so `--help` stays fast.
- **Never commit real data.** No real isolate IDs, patient or clinical data, assemblies, lab photos or
  workbooks. `data/`, `private/`, `outputs/` and `.env` are git-ignored. Keep them that way. Examples
  must stay fictional (`DEMO-*`) and synthetic.
- **No personal details or local paths.** No emails, phone numbers, institution-internal names or
  `C:\Users\...` paths in code, docs or file metadata.
- **Facts must be real.** Don't invent features, validation results, citations or performance claims in docs.
  `docs/verification.md` states what has and has not been validated. Keep it honest.
- **Keep it cross-platform.** Use `pathlib`, no shell-specific calls. External tools (BLAST+) are found on PATH.
- **Dependencies:** add new ones only when needed, to `pyproject.toml` and the matching `requirements*.txt`.

## Publishing

Remote (once created): `https://github.com/shadhin-ali002/amr-analysis-and-visualization`.

```powershell
git add -A
git commit -m "Describe the change"
git push
```

Commit, push, create the GitHub repo or change its visibility only when the owner asks.
