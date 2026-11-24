# Outputs

Each run creates `Run_<timestamp>/` under `--output` (default `outputs/`). Files marked *(genomic)* appear
only when `--wgs-base` is given and BLAST+ and the genomics extras are installed.

## `01_Mutation_Analysis/`

### `NG_AMR_MIC.xlsx`

| Sheet | Content |
|---|---|
| Isolate MICs | One row per isolate: every antibiotic × method MIC (µg/mL), cohort (clinical / QC) |
| Method details | Long format: MICs, S/I/R categories, McFarland and dates per isolate and antibiotic |
| ID reconciliation | Workbook IDs vs assembly files, and whether phenotype, sequence or both exist |
| QC | Results for the ATCC 49226 quality-control strain |
| Source cells | Sheet and cell each value came from, with its unit and parse status |
| Paired agreement | E-strip vs agar dilution: exact pairs, essential agreement (±1 doubling dilution), categorical agreement, mean log₂ bias, Wilcoxon signed-rank p |
| Recorded WGS agreement | Phenotype category vs the recorded WGS-predicted category |
| Method pairs | S/I/R agreement for every method pair, with exact McNemar p |
| AD comparator errors | Major, very major and minor errors of E-strip against agar dilution |

### Mutation detection *(genomic)*

| Path | Content |
|---|---|
| `Mutation_Tables/<gene>_protein_summary.csv` | Per isolate: identity, coverage, status (intact / disrupted / absent / partial), each variant with type, position and curated flag |
| `Mutation_Tables/<gene>_mutation_matrix.csv` | Isolate × variant presence (YES/NO) |
| `Mutation_Reports/By_Gene/`, `By_Isolate/` | PDF and HTML alignment reports with highlighted variants, plus Excel call tables |
| `NG_AMR_Mutation_MIC_Matrix.xlsx` | Mutation × isolate matrix per gene with MICs, call status and isolate metadata |
| `NG_AMR_MIC_MasterMatrix.xlsx` | Combined MIC and mutation master workbook |
| `WGS_calls.csv` | Status, identity and coverage of every locus in every isolate |
| `Extracted_Sequences/<isolate>/` | Called nucleotide and protein FASTA per gene |

### `Summary_Tables/` *(genomic)*

| Folder | Content |
|---|---|
| `00_NG_AMR_Master_Mutation_Matrix.xlsx` | Copy of the master matrix |
| `01_Mutation_Catalogue` | Every variant ranked by prevalence, with drug class and curated flag (CSV + XLSX) |
| `02_Mutation_Burden` | Mean variants per isolate and % isolates with ≥1 variant, per gene |
| `03_Threshold_Stratification` | Variants meeting each `--thresholds` prevalence cut-off |
| `04_Curated_vs_NonCurated` | Curated (catalogue) vs other variants per gene |
| `05_DrugClass_Summary` | Variants and events grouped by antibiotic class |
| `06_Key_Findings_and_Index` | `KEY_FINDINGS.md`: plain-text summary and folder index |
| `07_Typing_MLST_NGSTAR` | Template for typing results produced by external tools |
| `08_Phylogeny` | Checklist for external phylogeny |
| `09_Plasmid_MGE` | blaTEM / tetM detection summary |
| `10_Structural_penA_PBP2` | penA variants mapped to published PBP2 structural motifs |

## `05_Other_Outputs/`

| File | Content |
|---|---|
| `run_status.csv` | Every stage with `ok`, `skipped` or `failed` and a detail message |
| `NG_AMR_Run_Manifest.json` | Arguments, timings, Python and package version, SHA-256 of the workbook and of each source file |
| `analysis.log` | Full log of the run |

The command exits with status 1 if any stage failed.
