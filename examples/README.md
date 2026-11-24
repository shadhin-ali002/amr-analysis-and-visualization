# Examples

All inputs here are fictional.

| File | Content |
|---|---|
| `demo_mic.xlsx` | Five-sheet MIC workbook (`EST_MIC`, `AD_MIC`, `McFarland`, `Disc_Diffusion`, `WGS`) for isolates DEMO-A, DEMO-B and DEMO-C (penicillin G) |
| `demo_aliases.csv` | Maps the old ID `DEMO-OLD` to `DEMO-A` |

From the project folder:

```powershell
ng-amr --workbook examples/demo_mic.xlsx --aliases examples/demo_aliases.csv
```

For the genomic route, use the synthetic genome set in `tests/fixtures/genomic_demo`:

```powershell
ng-amr --workbook examples/demo_mic.xlsx --wgs-base tests/fixtures/genomic_demo
```
