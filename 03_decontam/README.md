# 03 — Contaminant removal (decontam)

Removes likely laboratory/reagent contaminants using the prevalence-based
method in the R `decontam` package, comparing true samples against
extraction/PCR negative controls.

## Input
- `phyloseq_MARES_full_e20.rds` (from `02_taxonomic_assignment/`)

## Script
- `run_decontam.R`

## Parameters
| Parameter | Value |
|---|---|
| Method | prevalence |
| Negative control samples | `N_N_N1`–`N_N_N5` |
| Threshold | 0.5 |

`decontam` version: `[TBD — see environment/]`

## Output
- `phyloseq_MARES_decontam_true_e20.rds` — contaminants removed, negative
  controls dropped
- `ASVtable_decontam_true_e20.tsv` — exported table, input for
  `04_tag_switch_correction/`
