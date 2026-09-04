# ══════════════════════════════════════════════════════════════════════════
# 02.3 — Build the ASV-level phyloseq object (OTU table + taxonomy +
#         metadata + representative sequences)
# ══════════════════════════════════════════════════════════════════════════
# INPUT:  ASVtable/ASV-frequency-table.tsv
#         taxonomy_edited_8ranks_e20.txt   (from 02_clean_taxonomy_path.sh)
#         sample_metadata_phyloseq.tsv
#         ASVseqs/rep-seq-ASV.fasta
# OUTPUT: phyloseq_MARES_full_e20.rds

library(phyloseq)
library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(Biostrings)

asv_tsv  <- "ASVtable/ASV-frequency-table.tsv"
tax_tsv  <- "taxonomy_edited_8ranks_e20.txt"
meta_tsv <- "sample_metadata_phyloseq.tsv"
fasta_asv <- if (file.exists("ASVseqs/rep-seq-ASV.fasta")) "ASVseqs/rep-seq-ASV.fasta" else "rep-seq-ASV.fasta"

# ---- 1) ASV table ----
hdr <- readLines(asv_tsv, n = 2)[2]
hdr <- sub("^#OTU ID\\s+", "", hdr)
samples <- str_split(hdr, "\\s+|\t+")[[1]]
samples <- samples[samples != ""]

tab <- read.delim(asv_tsv, sep = "\t", header = FALSE, skip = 2,
                   row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
colnames(tab) <- samples

otu_mat <- as.matrix(tab)
mode(otu_mat) <- "numeric"
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

# ---- 2) Taxonomy ----
tax_raw <- read_tsv(tax_tsv, show_col_types = FALSE)
stopifnot(all(c("OTU_ID","TaxonPath") %in% colnames(tax_raw)))

pick1 <- function(x) { x <- x[!is.na(x)]; if (length(x)==0) NA_character_ else x[1] }

tax_parsed <- tax_raw %>%
  mutate(TaxonPath = ifelse(is.na(TaxonPath) | TaxonPath=="", "d__unassigned;", TaxonPath),
         TaxonPath = str_replace_all(TaxonPath, "\\s+", " "),
         parts = str_split(TaxonPath, ";"),
         parts = lapply(parts, function(x) x[x!=""])) %>%
  rowwise() %>%
  mutate(
    Domain  = pick1(parts[str_starts(parts, "d__")]),
    Phylum  = pick1(parts[str_starts(parts, "p__")]),
    Class   = pick1(parts[str_starts(parts, "c__")]),
    Order   = pick1(parts[str_starts(parts, "o__")]),
    Family  = pick1(parts[str_starts(parts, "f__")]),
    Genus   = pick1(parts[str_starts(parts, "g__")]),
    Species = pick1(parts[str_starts(parts, "s__")])
  ) %>%
  ungroup() %>%
  select(OTU_ID, Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  mutate(across(Domain:Species, ~ str_remove(.x, "^[a-z]__")))

tax_mat <- as.matrix(tax_parsed %>% column_to_rownames("OTU_ID"))
TAX <- tax_table(tax_mat)

# ---- 3) Representative sequences ----
stopifnot(file.exists(fasta_asv))
REF <- readDNAStringSet(fasta_asv, format = "fasta")

# ---- 4) Sample metadata ----
meta <- read_tsv(meta_tsv, show_col_types = FALSE)
stopifnot("SampleID" %in% colnames(meta))
meta <- meta %>%
  mutate(SampleID = str_trim(SampleID),
         SampleID = str_replace_all(SampleID, "\\s+", "_"))
meta_df <- meta %>% column_to_rownames("SampleID")
SD <- sample_data(meta_df)

# ---- 5) Sync samples and taxa across all four components ----
common_samples <- intersect(sample_names(OTU), rownames(meta_df))
if (length(common_samples) == 0) {
  stop("No common samples — check that ASV table column names match metadata SampleID.")
}

OTU2 <- prune_samples(common_samples, OTU)
SD2  <- sample_data(meta_df[common_samples, , drop = FALSE])

common_taxa <- intersect(taxa_names(OTU2), taxa_names(TAX))
OTU2 <- prune_taxa(common_taxa, OTU2)
TAX2 <- tax_table(TAX)[common_taxa, ]

common_taxa2 <- intersect(common_taxa, names(REF))
OTU2 <- prune_taxa(common_taxa2, OTU2)
TAX2 <- tax_table(TAX2)[common_taxa2, ]
REF2 <- REF[common_taxa2]

# ---- 6) Build phyloseq ----
ps <- phyloseq(OTU2, SD2, TAX2, refseq(REF2))
saveRDS(ps, "phyloseq_MARES_full_e20.rds")
ps
