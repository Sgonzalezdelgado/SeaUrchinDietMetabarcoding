# ══════════════════════════════════════════════════════════════════════════
# 04.2 — Rebuild phyloseq object after tag-switch correction
# ══════════════════════════════════════════════════════════════════════════
# INPUT:  phyloseq_MARES_decontam_true_e20.rds  (from 03_decontam/)
#         ASVtable_tsc_e20.tsv                  (from run_owi_renormalize.sh)
# OUTPUT: phyloseq_MARES_decontam_tsc_e20.rds
#         (this is the file used as input everywhere downstream)

library(phyloseq)
library(Biostrings)

ps_true <- readRDS("phyloseq_MARES_decontam_true_e20.rds")

tab <- read.delim(
  "ASVtable_tsc_e20.tsv",
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
)

# owi_renormalize adds extra non-sample columns (e.g. "Total") — drop them
tab <- tab[, colnames(tab) %in% sample_names(ps_true), drop = FALSE]
tab <- tab[rowSums(tab) > 0, , drop = FALSE]

OTU_tsc <- otu_table(as.matrix(tab), taxa_are_rows = TRUE)

common_samp <- intersect(sample_names(ps_true), colnames(tab))
common_taxa <- intersect(taxa_names(ps_true), rownames(tab))

ps_tsc <- phyloseq(
  prune_taxa(common_taxa, OTU_tsc),
  sample_data(prune_samples(common_samp, ps_true)),
  tax_table(ps_true)[common_taxa, ],
  refseq(ps_true)[common_taxa]
)

saveRDS(ps_tsc, "phyloseq_MARES_decontam_tsc_e20.rds")
ps_tsc
