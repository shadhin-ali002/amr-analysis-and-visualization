# AGENTS.md

Guidance for AI coding agents working in this repository.

## What this is

`amr-analysis-and-visualization` is a Python command-line tool for antimicrobial-resistance (AMR) work on
*Neisseria gonorrhoeae* isolates, published by Md. Shadhin Ali (`shadhin-ali002`). It:

- reads an MIC workbook (sheets `EST_MIC`, `AD_MIC`, `McFarland`, `Disc_Diffusion`, `WGS`) and writes MIC
  tables, figures, E-strip vs agar-dilution comparisons and a run-summary presentation;
- optionally, with existing genome assemblies and reference files, calls mutations in resistance genes with
  BLAST+ and writes grouped reports by gene and by isolate (PDF, HTML, Excel).

It is one of five sibling repos (`e-strip-plate-analysis`, `agar-dilution-plate-analysis`,
`phylogeny-analysis`, `virulence-analysis`). Each must run on its own. The owner's profile README and
website live in `../shadhin-ali002` and `../shadhin-ali002.github.io`.

## Files

| Path | Purpose |
| --- | --- |
| `run_amr.py` | CLI entry point: argument parsing and the run sequence |
| `amr/inputs.py` | ID normalisation (`canonical`, aliases, QC strain ATCC 49226), `Audit`, assembly discovery, workbook constants |
| `amr/engine.py` | Most of the logic (~8.7k lines): workbook reading, BLAST and mutation calling, charts, reports, `run_self_tests()` |
| `amr/mutation_reports.py` | Grouped by-gene / by-isolate mutation reports |
| `amr/visualization.py` | Chart helpers |
| `amr/phenotype.py`, `amr/mutations.py`, `amr/run_summary.py` | Thin stage wrappers called by `run_amr.py` |
| `examples/` | Fictional demo workbook (`DEMO-A/B/C`, penicillin G) and alias CSV |
| `tests/test_example.py` | End-to-end test of the example run |
| `docs/` | Input formats, verification notes |
| `LICENSE` | MIT. Don't change it, or add code under an incompatible licence, without the owner's request |
| `requirements.txt` / `requirements-genomics.txt` | Core deps / extra deps for the genomic route (Biopython, parasail) |

A restructure into a `src/ng_amr/` package with sub-packages is planned. When it lands, update this table.

## Setup and checks

Windows PowerShell, Python 3.11+:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-genomics.txt
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
.\.venv\Scripts\python.exe run_amr.py --self-test
.\.venv\Scripts\python.exe run_amr.py --workbook examples/demo_mic.xlsx --output $env:TEMP\amr-out
```

Expected: the unit test passes, the self-test reports `ALL SELF-TESTS PASSED (33 checks)` (23 without the
genomics extras), and the example run exits 0 with `run_status.csv` showing no `failed` rows.
The genomic route also needs NCBI BLAST+ (`blastn`, `makeblastdb`) on PATH.

Run all three checks after every change. Write example output outside the repo or to `outputs/` (git-ignored).

## Rules

- **Do not change results.** This produces scientific numbers. Refactors must keep every output file name,
  folder name, column header and value identical. Compare the example output before and after a change.
  Changes to calling logic, thresholds, breakpoints or the curated catalogue need the owner's explicit request.
- **Move code, don't rewrite it.** When restructuring, cut and paste functions unchanged. Make one move at a
  time and re-run the checks.
- **Never commit real data.** No real isolate IDs, patient or clinical data, assemblies, lab photos or
  workbooks. `data/`, `private/`, `outputs/` and `.env` are git-ignored. Keep them that way. Examples
  must stay fictional (`DEMO-*`) and synthetic.
- **No personal details or local paths.** No emails, phone numbers, institution-internal names or
  `C:\Users\...` paths in code, docs or file metadata.
- **Facts must be real.** Don't invent features, validation results, citations or performance claims in docs.
  `docs/verification.md` states what has and has not been validated. Keep it honest.
- **Keep it cross-platform.** Use `pathlib`, no shell-specific calls. External tools (BLAST+) are found on PATH.
- **Dependencies:** add new ones only when needed, and to both `requirements*.txt` (and `pyproject.toml` once it exists).

## Publishing

Remote (once created): `https://github.com/shadhin-ali002/amr-analysis-and-visualization`.

```powershell
git add -A
git commit -m "Describe the change"
git push
```

Commit, push, create the GitHub repo or change its visibility only when the owner asks.
