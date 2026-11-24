# Input formats

## MIC workbook (`--workbook`)

An `.xlsx` file with the five sheets of `examples/demo_mic.xlsx`: `EST_MIC`, `AD_MIC`, `McFarland`,
`Disc_Diffusion` and `WGS`. Keep the isolate-ID headers and the unit in MIC headers (for example
`Penicillin G MIC (ug/mL)`). E-strip and agar-dilution values go in separate fields. MIC values may carry a
bound (`<=0.016`, `>256`). If cells are formula-derived, save the workbook so that the calculated values are
stored.

## Genome data (`--wgs-base`)

```
<wgs-base>/
├── Assembly/
│   ├── NT/Plane_Fasta/<isolate>.fasta   assembly contigs, one file per isolate (required)
│   ├── NT/<isolate>.fasta               annotated nucleotide FASTA (optional)
│   └── AA/                              original protein FASTA for comparison (optional)
└── Refs/
    ├── Genes_RefSeq.fasta               one reference sequence per locus, FASTA ID = locus name
    └── Curated_Variants.json            {"gyrA": ["S91F", ...], "penA": [...], ...}
```

Locus names must match the panel in `src/ng_amr/core/constants.py` (`DEFAULT_GENES`). File names become
isolate IDs: a leading number and underscore are dropped, and text is upper-cased. `tests/fixtures/genomic_demo`
is a complete synthetic example.

## Optional CSV files

| Option | Header | Purpose |
|---|---|---|
| `--aliases` | `alias,canonical_id` | Map each historical ID directly to its final ID |
| `--id-map` | `filename,isolate_id` | Give the isolate ID for an assembly file name |

Keep real inputs in the git-ignored `data/` or `private/` folders.
