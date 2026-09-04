# 04 — Tag-switching correction

Correction of tag-jump/tag-switching artefacts with `owi_renormalize`
(Wangensteen, Project Metabarpark — https://github.com/metabarpark/R_scripts_metabarpark).

## Input
- `phyloseq_MARES_decontam_tsc_e20.rds` (from step 01)

## Script
- `owi_renormalize.R`

## Command used
```bash
Rscript owi_renormalize.R -c 0.97 -s <start_col> -e <end_col>
```
`-s`/`-e` are the start/end **column indices** of the sample columns in the
exported table — not method parameters. Only `-c` (cumulative frequency
cutoff = 0.97) is scientifically meaningful and needs to match the Methods
text.

## Output
- Renormalized OTU/ASV table, feeding into step 05
