# Input formats

| Input | Meaning |
|---|---|
| MIC workbook | Use the five-sheet example in the examples folder. Preserve the isolate-ID headers and the unit in MIC headers. Stored E-strip and agar-dilution observations have distinct fields. |
| Private ID aliases | CSV header `alias,canonical_id`. Each alias maps directly to one final identifier. See the fictional example in examples/demo_aliases.csv. |
| Assembly filename map | CSV header `filename,isolate_id`. Use with `--id-map` when supplying assemblies. |

Local private inputs, aliases and outputs can go under the ignored private/, data/ and outputs/ directories.
