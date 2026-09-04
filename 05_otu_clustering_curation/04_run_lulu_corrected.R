# ══════════════════════════════════════════════════════════════════════════
# 03.4 — Run corrected LULU curation
# ══════════════════════════════════════════════════════════════════════════
# INPUT:  phyloseq_OTU97_decontam_tsc_e20.rds  (from 02_build_otu97_phyloseq.R)
#         match_list_OTU97_3col_e20.txt        (from 03_selfmatch_vsearch.sh)
#         lulu2_function.R                     (patched function, this folder)
# OUTPUT: phyloseq_OTU97_LULU_e20.rds  (curated: 30,386 taxa x 170 samples)

library(phyloseq)

source("lulu2_function.R")   # defines lulu2() in this session

ps <- readRDS("phyloseq_OTU97_decontam_tsc_e20.rds")

otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu <- t(otu)
mode(otu) <- "numeric"
otu <- otu[rowSums(otu) > 0, , drop = FALSE]
otu_df <- as.data.frame(otu, check.names = FALSE)

match_raw <- read.delim(
  "match_list_OTU97_3col_e20.txt",
  header = FALSE, sep = "\t",
  stringsAsFactors = FALSE
)
colnames(match_raw) <- c("query", "target", "identity")
match_raw$identity <- as.numeric(match_raw$identity)
if (max(match_raw$identity, na.rm = TRUE) > 1) {
  match_raw$identity <- match_raw$identity / 100
}

match_df <- data.frame(
  OTU = as.character(match_raw$query),
  Match = as.character(match_raw$target),
  Similarity = as.numeric(match_raw$identity),
  stringsAsFactors = FALSE
)
match_df <- match_df[!is.na(match_df$Similarity), ]
match_df <- match_df[match_df$OTU != match_df$Match, ]
keep <- match_df$OTU %in% rownames(otu_df) & match_df$Match %in% rownames(otu_df)
match_df <- match_df[keep, , drop = FALSE]

# minimum_match = 0.84 (not the default 84) because Similarity is scaled 0-1
cur <- lulu2(otu_df, match_df,
             minimum_ratio_type = "min",
             minimum_relative_cooccurence = 0.5,
             minimum_match = 0.84)

cat("Curated OTUs:", cur$curated_count, "\n")
cat("Discarded (merged) OTUs:", cur$discarded_count, "\n")

# ---- rebuild phyloseq with the curated table ----
otu_cur <- cur$curated_table
OTUcur <- otu_table(as.matrix(otu_cur), taxa_are_rows = TRUE)

common_taxa <- intersect(rownames(otu_cur), taxa_names(ps))

ps_lulu <- phyloseq(
  prune_taxa(common_taxa, OTUcur),
  sample_data(ps),
  tax_table(ps)[common_taxa, ],
  refseq(ps)[common_taxa]
)

saveRDS(ps_lulu, "phyloseq_OTU97_LULU_e20.rds")
ps_lulu
