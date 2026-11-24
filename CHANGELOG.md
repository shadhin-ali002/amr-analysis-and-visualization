# Changelog

## 1.0.0 — 2026-11-25

First public release.

- Code organised as the installable package `ng_amr` (`src/` layout) with sub-packages `core`, `io`,
  `phenotype`, `genomics` and `reports`. Installs an `ng-amr` command; `python run_amr.py` still works.
- Outputs are tables (XLSX/CSV), reports (PDF/HTML) and text only. All PNG figure, PowerPoint and Word
  gallery generation was removed, along with the matplotlib, Pillow, python-pptx and python-docx
  dependencies. A full genomic run of the demo now takes seconds instead of minutes.
- The publication tables folder was renamed to `Summary_Tables/` and its sub-folders renumbered.
- The run manifest records the package version and checksums of every source module.
- New tests: CLI argument validation, identifier and workbook parsing, built-in checks under `unittest`, and
  an end-to-end BLAST+ test on a synthetic genome fixture with planted variants.
- Added `pyproject.toml`, licence (MIT), citation file and a GitHub Actions test workflow.
