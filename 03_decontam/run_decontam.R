# ══════════════════════════════════════════════════════════════════════════
# 03 — Contaminant removal (decontam, prevalence method)
# ══════════════════════════════════════════════════════════════════════════
# INPUT:  phyloseq_MARES_full_e20.rds
# OUTPUT: phyloseq_MARES_decontam_true_e20.rds
#         ASVtable_decontam_true_e20.tsv

library(phyloseq)
library(decontam)

ps <- readRDS("phyloseq_MARES_full_e20.rds")

# Adjust to the exact negative control sample IDs as they appear in
# sample_names(ps)
neg_ids <- c("N_N_N1", "N_N_N2", "N_N_N3", "N_N_N4", "N_N_N5")

sample_data(ps)$is.neg <- sample_names(ps) %in% neg_ids
cat("Negative controls found:", sum(sample_data(ps)$is.neg), "\n")

contam.prev <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.5)

ps_nocon <- prune_taxa(!contam.prev$contaminant, ps)
ps_true  <- subset_samples(ps_nocon, is.neg == FALSE)

saveRDS(ps_true, "phyloseq_MARES_decontam_true_e20.rds")

# Export table for the tag-switch correction step
otu <- as(otu_table(ps_true), "matrix")
if (!taxa_are_rows(ps_true)) otu <- t(otu)

write.table(otu, file = "ASVtable_decontam_true_e20.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

cat("OK: phyloseq_MARES_decontam_true_e20.rds y ASVtable_decontam_true_e20.tsv\n")
