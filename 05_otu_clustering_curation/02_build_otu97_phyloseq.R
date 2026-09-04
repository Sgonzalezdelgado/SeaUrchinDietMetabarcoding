# ══════════════════════════════════════════════════════════════════════════
# 03.2 — Build the OTU97 phyloseq object
# ══════════════════════════════════════════════════════════════════════════
# INPUT:  phyloseq_MARES_decontam_tsc_e20.rds (ASV-level, from step 02)
#         OTU97_table_e20.tsv, ASV_to_OTU97_e20.tsv, OTU97_centroids_e20.fasta
#         (from 01_cluster_vsearch.sh)
# OUTPUT: phyloseq_OTU97_decontam_tsc_e20.rds  (32,993 taxa x 170 samples)

library(phyloseq)
library(dplyr)
library(readr)
library(tibble)
library(Biostrings)

ps_asv <- readRDS("phyloseq_MARES_decontam_tsc_e20.rds")

asv_tab <- read.delim("OTU97/ASV_table_tsc_e20.tsv", sep = "\t",
                       header = TRUE, row.names = 1, check.names = FALSE)
map <- read_tsv("OTU97/ASV_to_OTU97_e20.tsv",
                 col_names = c("ASV", "OTU97"), show_col_types = FALSE)

asv_tab$ASV <- rownames(asv_tab)
df <- asv_tab %>% left_join(map, by = "ASV")
df$OTU97[is.na(df$OTU97)] <- df$ASV[is.na(df$OTU97)]

otu97_tab <- df %>%
  select(-ASV) %>%
  group_by(OTU97) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  column_to_rownames("OTU97")

write.table(otu97_tab, "OTU97/OTU97_table_e20.tsv", sep = "\t",
            quote = FALSE, col.names = NA)

OTU97 <- otu_table(as.matrix(otu97_tab), taxa_are_rows = TRUE)

centroids <- readDNAStringSet("OTU97/OTU97_centroids_e20.fasta", format = "fasta")

tax_asv <- as(tax_table(ps_asv), "matrix")
otu_ids <- taxa_names(OTU97)
keep_tax <- intersect(otu_ids, rownames(tax_asv))
tax_otu <- tax_asv[keep_tax, , drop = FALSE]

missing <- setdiff(otu_ids, rownames(tax_otu))
if (length(missing) > 0) {
  fill <- matrix(NA_character_, nrow = length(missing), ncol = ncol(tax_asv))
  colnames(fill) <- colnames(tax_asv)
  rownames(fill) <- missing
  tax_otu <- rbind(tax_otu, fill)
}
tax_otu <- tax_otu[otu_ids, , drop = FALSE]
TAX97 <- tax_table(tax_otu)

SD <- sample_data(ps_asv)

common_otu <- intersect(otu_ids, names(centroids))
OTU97 <- prune_taxa(common_otu, OTU97)
TAX97 <- tax_table(TAX97)[common_otu, ]
centroids2 <- centroids[common_otu]

ps_otu97 <- phyloseq(OTU97, SD, TAX97, refseq(centroids2))

saveRDS(ps_otu97, "phyloseq_OTU97_decontam_tsc_e20.rds")
ps_otu97
