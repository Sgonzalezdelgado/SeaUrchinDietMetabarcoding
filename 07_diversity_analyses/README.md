# 07 — Diversity and statistical analyses

Three objective-specific scripts, all starting from the final curated
phyloseq object produced in `06_taxonomic_curation/`

| Script | Objective |
|---|---|
| `O1_KKvsE_final.R` | Validation of faecal vs. stomach content sampling (Blanes) |
| `O2_VentvsControl_final.R` | Diet composition along the pH gradient (La Palma vent system) |
| `O3_alixulavsPlividus_final.R` | Interspecific dietary overlap (*A. lixula* vs *P. lividus*) |

Each core script covers: alpha diversity (observed richness, LMM),
beta diversity (Jaccard dissimilarity, PERMANOVA via `adonis2`), and,
where relevant, indicator taxon analysis (`indicspecies`) and dietary
overlap. Venn diagrams in O1 are hand-built with `ggplot2`
(`geom_polygon`)

## Note
Alpha diversity: test pH condition while accounting for locality as a random effect.
Beta diversity: test locality only, as condition-level permutations are not robust with 
only five independent localities (two vent and three ambient).

## Input
- Final curated phyloseq object from `06_taxonomic_curation/`
  (3,686 OTUs / 838 species-level taxa)

## Output
- Figures 2–5, Tables S2–S10 (all confirmed against the corrected LULU
  curation as of the last sync)
