# 06 — Taxonomic assignment and dietary curation

## Steps
1. Taxonomic assignment with MEGAN6 (min score 98, max expected 1×10⁻⁸,
   min identity 80%, top 10%)
2. Restriction to Eukaryota.
3. **Rule-based curation** (`precuration.R`, 29 documented steps): removes
   terrestrial, freshwater, parasitic, out-of-range, and lab-contaminant
   taxa by phylum/genus/species-specific rules. Not manual clicking —
   fully scripted, so re-running it on a new input table is mechanical.
   Host reads (*A. lixula* / *P. lividus*) are **not** removed at this
   step; they're stripped later, per-objective, in `07_diversity_analyses/`
4. Species-level collapse (`tax_glom`) of OTUs sharing the same taxonomic
   assignment.

## Input
- `phyloseq_OTU97_LULU_e20.rds` (corrected, 30,386 OTUs — from step 05)

## Output
- Final phyloseq object used for all diversity analyses (Objectives 1–3):
  3,686 OTUs / 838 species-level taxa / 141 KK+E samples (La Palma) after
  further per-objective filtering in `07_diversity_analyses/`.