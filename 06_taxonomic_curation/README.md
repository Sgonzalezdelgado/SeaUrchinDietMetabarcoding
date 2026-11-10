# 06 — Taxonomic assignment, dietary curation and seawater decontamination

## Steps
| Script | What it does |
|---|---|
| `precuration.R` | Eukaryote restriction and rule-based dietary curation (documented steps). Writes the removal log and the supplementary OTU table and collapses OTUs to species level. |
| `water_microdecon.R` | **Main seawater decontamination.** Uses the seawater samples (tissue `W`) as blanks in microDecon (McKnight et al. 2019, *Environmental DNA* 1:14–25). |
| `water_overlap_filter.R` | Sensitivity analysis only (`MODE = "presence"` removes every OTU found in the paired water). Not used for the main results. |

## Input
- `phyloseq_OTU97_LULU_e20.rds` (from step 05; includes faeces, stomachs and seawater)

## Output
| File | Used by |
|---|---|
| `phyloseq_OTU97_LULU_e20_EUK_KKE_curated.rds` | `water_microdecon.R` |
| `phyloseq_OTU97_LULU_e20_EUK_KKE_curated_speciesGlom.rds` | `O1_KKvsE.R` (Blanes) |
| `..._curated_waterFiltered_microdecon.rds` | — (OTU level, archived) |
| `..._curated_waterFiltered_microdecon_speciesGlom.rds` | `O2_VentvsAmbient.R`, `O3_alixulavsPlividus.R` |
| `microdecon_*.csv`, `S1_*.csv` | Supplementary Tables S1 and S11 |
