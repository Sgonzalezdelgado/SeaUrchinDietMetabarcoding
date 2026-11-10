# 07 — Diversity and statistical analyses

Three objective-specific scripts, all starting from the final curated
phyloseq object produced in `06_taxonomic_curation/`

| Script | Objective |
|---|---|
| `O1_KKvsE_final.R` | Validation of faecal vs. stomach content sampling (Blanes) |
| `O2_VentvsControl_final.R` | Diet composition along the pH gradient (La Palma vent system) |
| `O3_alixulavsPlividus_final.R` | Interspecific dietary overlap (*A. lixula* vs *P. lividus*) |
| `plot_indval_robust.R` | Figure 4 from the robust IndVal results |

Each core script covers: alpha diversity (observed richness, LMM),
beta diversity (Jaccard dissimilarity, PERMANOVA via `adonis2`), and,
where relevant, indicator taxon analysis (`indicspecies`) and dietary
overlap. Venn diagrams in O1 are hand-built with `ggplot2`
(`geom_polygon`)

## Input
- Final curated phyloseq object from `06_taxonomic_curation/`
  (3,686 OTUs / 838 species-level taxa)

## Output
- Objective 1: `OBJ1_*` files (Figure 2, Figure S1, Tables S2–S3).
- Objective 2: `Objetivo2_water_microdecon/` (Figure 3, Figure 4, Tables S4–S7).
- Objective 3: `Objetivo3_water_microdecon/` (Figure 5, Figure S3, Tables S8–S10).
