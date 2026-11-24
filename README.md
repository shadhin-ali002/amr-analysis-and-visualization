# AMR analysis and visualization

Summarize observed MICs and paired methods, and produce grouped mutation reports by gene and isolate when existing assemblies and reference files are available. Plain workbook mode reads the five sheets shown in the example: EST_MIC, AD_MIC, McFarland, Disc_Diffusion and WGS. Units and isolate headers matter. Supply a workbook with saved numerical values for formula-derived cells.

The example contains three fictional samples. Inspect the “Isolate MICs” sheet of `01_Mutation_Analysis/NG_AMR_MIC.xlsx`. The example AD column for penicillin G contains 0.25, 0.5 and 1 µg/mL. The input categories are fictional.

For the optional genomic workflow install `requirements-genomics.txt` and place NCBI BLAST+ executables `blastn` and `makeblastdb` on PATH. Use `--wgs-base` with a prepared dataset organized as `Assembly/NT/Plane_Fasta/` for assembly FASTAs, `Assembly/NT/` for annotated inputs, and optionally `Assembly/AA/`. The reference directory contains `Refs/Genes_RefSeq.fasta` and `Refs/Curated_Variants.json`. Select an existing catalogue of another name with `--catalog`. Supply these data separately. Grouped PDF, HTML and Excel outputs reside under `Mutation_Reports/By_Gene/` and `Mutation_Reports/By_Isolate/` inside the mutation output folder. Use `--id-map` for an explicit filename,isolate_id CSV. Run `python run_amr.py --self-test` for the inherited regression checks.

Generated figures, workbooks and grouped reports are under `01_Mutation_Analysis/`. The run log, manifest, stage status and summary presentation are under `05_Other_Outputs/`.

## First run on Windows

Open PowerShell in this project folder. Use Python 3.11 or newer; choose the corresponding installed interpreter in the first command.

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe run_amr.py --workbook examples/demo_mic.xlsx
```

Use `--output` for a custom output directory. The default is `outputs/` in the working directory. Run `python run_amr.py --help` for all flags. Run `python -m unittest discover -s tests -v` with the project’s environment to check the supplied example. On macOS or Linux create the environment with `python3 -m venv .venv` and use `.venv/bin/python`.

[Input formats](docs/data-formats.md) · [Verification](docs/verification.md) · [Upload to GitHub](docs/github-upload.md)
