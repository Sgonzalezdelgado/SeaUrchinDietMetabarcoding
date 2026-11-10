# ══════════════════════════════════════════════════════════════════════════════
# WATER DECONTAMINATION WITH microDecon — run after precuration.R
# ══════════════════════════════════════════════════════════════════════════════
# Faecal pellets were collected by pipette from tanks filled with seawater from
# each locality. Here the seawater samples of each locality are used as
# "blanks" in microDecon (McKnight et al. 2019, Environmental DNA 1:14-25),
# which subtracts water-derived reads OTU by OTU, proportionally to the water
# composition, instead of removing whole OTUs.
#
# Design choices (see microDecon user's guide, sections 1.4.4-1.4.5):
#   - One run per batch (locality x campaign), each with its own water as blanks,
#     because the contamination source (water) differs among localities.
#   - Within each batch, host species are the a priori groups (numb.ind) used
#     by the residual-contamination filter (thresh, prop.thresh).
#   - Host and Echinodermata OTUs removed before decontamination (as in O2/O3).
#   - Blanes has no water samples and is left unchanged.
#   - microDecon defaults otherwise (runs = 2, thresh = 0.7, prop.thresh = 5e-5).
#
# INPUTS:
#   phyloseq_OTU97_LULU_e20.rds                  (all samples, incl. water "W")
#   phyloseq_OTU97_LULU_e20_EUK_KKE_curated.rds  (OTU level, precuration Step 27)
# OUTPUTS:
#   microdecon_water_heterogeneity.csv   Bray-Curtis between water replicates
#   microdecon_per_sample.csv            reads/OTUs removed per faecal sample
#   microdecon_removed_by_taxon.csv      reads removed per OTU, with taxonomy
#   microdecon_water_taxa_prevalence.csv prevalence in faeces of dominant water OTUs
#   phyloseq_OTU97_LULU_e20_EUK_KKE_curated_waterFiltered_microdecon.rds
#   phyloseq_OTU97_LULU_e20_EUK_KKE_curated_waterFiltered_microdecon_speciesGlom.rds
# ══════════════════════════════════════════════════════════════════════════════

# install once: remotes::install_github("donaldtmcknight/microDecon")
library(microDecon)
library(phyloseq)
library(vegan)
library(dplyr)

setwd(".")
set.seed(42)

# ── PARAMETERS ────────────────────────────────────────────────────────────────
# Faecal locality -> water locality. NA = no water, left unchanged.
WATER_LOCALITY_MAP <- c("Transition Faro" = "Los Porretos",
                        "Blanes"          = NA)
SUFFIX <- "waterFiltered_microdecon"

# ── LOAD ──────────────────────────────────────────────────────────────────────
ps_all <- readRDS("phyloseq_OTU97_LULU_e20.rds")
ps_cur <- readRDS("phyloseq_OTU97_LULU_e20_EUK_KKE_curated.rds")

as_mat <- function(ps) { m <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) m <- t(m); m }

ps_w  <- subset_samples(ps_all, tissue == "W")
sd_w  <- as(sample_data(ps_w), "data.frame")
sd_c  <- as(sample_data(ps_cur), "data.frame")
otu_c <- as_mat(ps_cur)
otu_w <- as_mat(ps_w)
stopifnot(all(rownames(otu_c) %in% rownames(otu_w)))
otu_w <- otu_w[rownames(otu_c), , drop = FALSE]

cat("Water samples per locality:\n"); print(table(sd_w$locality))

# Remove host + Echinodermata (same rule as O2/O3)
tax_c    <- as.data.frame(tax_table(ps_cur))
host_idx <- grepl("Arbacia|Paracentrotus", tax_c$Species, ignore.case = TRUE) |
            (!is.na(tax_c$Phylum) & tax_c$Phylum == "Echinodermata")
otu_c <- otu_c[!host_idx, , drop = FALSE]
otu_w <- otu_w[!host_idx, , drop = FALSE]
tax_c <- tax_c[!host_idx, , drop = FALSE]
cat("Host/Echinodermata OTUs excluded:", sum(host_idx), "\n")

# ── 1. HETEROGENEITY AMONG WATER REPLICATES (guide section 1.4.5) ─────────────
het <- bind_rows(lapply(unique(sd_w$locality), function(loc) {
  ids <- rownames(sd_w)[sd_w$locality == loc]
  if (length(ids) < 2) return(data.frame(locality = loc, pair = NA, BrayCurtis = NA))
  m  <- t(otu_w[, ids, drop = FALSE]); m <- m / pmax(rowSums(m), 1)
  bc <- as.matrix(vegdist(m, method = "bray"))
  pr <- combn(ids, 2)
  data.frame(locality = loc, pair = paste(pr[1, ], pr[2, ], sep = " vs "),
             BrayCurtis = round(bc[t(pr)], 3))
}))
write.csv(het, "microdecon_water_heterogeneity.csv", row.names = FALSE)
cat("\nWater heterogeneity (Bray-Curtis on proportions):\n"); print(het)

# ── 2. microDecon PER LOCALITY ────────────────────────────────────────────────
kk_ids   <- rownames(sd_c)[sd_c$tissue == "KK"]
kk_loc   <- as.character(sd_c[kk_ids, "locality"])
kk_wloc  <- ifelse(kk_loc %in% names(WATER_LOCALITY_MAP), WATER_LOCALITY_MAP[kk_loc], kk_loc)

# Batches = locality x sampling campaign. Water from the two campaigns differs
# strongly at La Bajita and El Faro (Bray-Curtis 0.71-0.94 between campaigns),
# so each campaign's faeces are decontaminated with that campaign's water only.
# Campaign rule: faecal replicate numbers 1-30 = first campaign,
# 50+ = second campaign. Water CB_ and VF_ = first campaign; all other water =
# second campaign.
OLD_WATER <- "^(CB|VF)_"
w_camp  <- ifelse(grepl(OLD_WATER, rownames(sd_w)), "first", "second")
w_batch <- setNames(paste(sd_w$locality, w_camp), rownames(sd_w))

rep_num <- as.numeric(sub("^[^_]+_([0-9]+).*$", "\\1", kk_ids))   # PV_4A -> 4
stopifnot(!any(is.na(rep_num)))
kk_camp  <- ifelse(rep_num < 50, "first", "second")
kk_batch <- ifelse(is.na(kk_wloc), NA, paste(kk_wloc, kk_camp))
cat("\nFaecal samples per batch:\n"); print(table(kk_batch, useNA = "ifany"))
cat("Water samples per batch:\n");   print(table(w_batch))
otu_filt <- otu_c
removed_list <- list()

for (wl in sort(unique(na.omit(kk_batch)))) {
  w_ids <- names(w_batch)[w_batch == wl]
  if (length(w_ids) == 0) { warning("No water for ", wl); next }

  s_ids <- kk_ids[!is.na(kk_batch) & kk_batch == wl]
  sp    <- as.character(sd_c[s_ids, "host_species"])
  s_ids <- s_ids[order(sp)]                     # samples ordered by group
  n_ind <- as.vector(table(sort(sp)))

  mat  <- cbind(otu_w[, w_ids, drop = FALSE], otu_c[, s_ids, drop = FALSE])
  keep <- rowSums(mat) > 0
  dat  <- data.frame(OTU = rownames(mat)[keep], mat[keep, , drop = FALSE],
                     check.names = FALSE, stringsAsFactors = FALSE)

  cat("\n== ", wl, ": ", length(w_ids), " water | ",
      paste(names(table(sort(sp))), n_ind, collapse = ", "), " ==\n", sep = "")

  res <- decon(data = dat, numb.blanks = length(w_ids), numb.ind = n_ind, taxa = FALSE)

  # rebuild full matrix for these samples (OTUs totally removed are absent)
  dt  <- res$decon.table
  dt[is.na(dt)] <- 0
  stopifnot(all(s_ids %in% colnames(dt)))
  new <- matrix(0, nrow = nrow(otu_c), ncol = length(s_ids),
                dimnames = list(rownames(otu_c), s_ids))
  new[as.character(dt$OTU), s_ids] <- as.matrix(dt[, s_ids])
  otu_filt[, s_ids] <- new

  rr <- res$reads.removed
  rr$locality <- wl
  removed_list[[wl]] <- rr
}

# ── 3. SUMMARIES ──────────────────────────────────────────────────────────────
per_sample <- data.frame(
  sample        = kk_ids,
  locality      = kk_loc,
  batch         = kk_batch,
  host_species  = sd_c[kk_ids, "host_species"],
  reads_before  = colSums(otu_c[, kk_ids]),
  reads_after   = colSums(otu_filt[, kk_ids]),
  OTUs_before   = colSums(otu_c[, kk_ids] > 0),
  OTUs_after    = colSums(otu_filt[, kk_ids] > 0)
) %>% mutate(pct_reads_removed = round(100 * (1 - reads_after / pmax(reads_before, 1)), 2),
             pct_OTUs_removed  = round(100 * (1 - OTUs_after / pmax(OTUs_before, 1)), 2))
write.csv(per_sample, "microdecon_per_sample.csv", row.names = FALSE)

cat("\nRemoved per locality and species (median %):\n")
print(per_sample %>% group_by(locality, host_species) %>%
        summarise(n = n(), med_pct_reads = median(pct_reads_removed),
                  max_pct_reads = max(pct_reads_removed),
                  med_pct_OTUs = median(pct_OTUs_removed), .groups = "drop"), n = 40)

diff_mat  <- otu_c[, kk_ids] - otu_filt[, kk_ids]
by_taxon <- tax_c[, c("Phylum", "Class", "Family", "Genus", "Species")] %>%
  mutate(OTU = rownames(otu_c),
         reads_KK        = rowSums(otu_c[, kk_ids]),
         reads_removed   = rowSums(diff_mat),
         pct_removed     = round(100 * reads_removed / pmax(reads_KK, 1), 1),
         n_samples_zeroed = rowSums(otu_c[, kk_ids] > 0 & otu_filt[, kk_ids] == 0)) %>%
  filter(reads_removed > 0) %>% arrange(desc(reads_removed))
write.csv(by_taxon, "microdecon_removed_by_taxon.csv", row.names = FALSE)
cat("\nReads removed by phylum:\n")
print(by_taxon %>% group_by(Phylum) %>%
        summarise(n_OTUs = n(), reads_removed = sum(reads_removed),
                  reads_KK = sum(reads_KK), .groups = "drop") %>%
        mutate(pct = round(100 * reads_removed / reads_KK, 1)) %>%
        arrange(desc(reads_removed)), n = 30)

# ── 4. PREVALENCE IN FAECES OF THE DOMINANT WATER OTUs ────────────────────────
# If faecal DNA came from the water, the most abundant water taxa should appear
# in most faecal samples of the same locality.
prev <- bind_rows(lapply(unique(na.omit(kk_wloc)), function(wl) {
  w_ids <- rownames(sd_w)[sd_w$locality == wl]
  s_ids <- kk_ids[!is.na(kk_wloc) & kk_wloc == wl]
  if (length(w_ids) == 0) return(NULL)
  wp  <- rowMeans(sweep(otu_w[, w_ids, drop = FALSE], 2,
                        pmax(colSums(otu_w[, w_ids, drop = FALSE]), 1), "/"))
  top <- names(sort(wp, decreasing = TRUE))[1:min(20, sum(wp > 0))]
  data.frame(locality = wl, OTU = top,
             pct_in_water   = round(100 * wp[top], 2),
             prev_faeces    = round(100 * rowMeans(otu_c[top, s_ids, drop = FALSE] > 0), 1),
             tax_c[top, c("Phylum", "Class", "Genus", "Species")])
}))
write.csv(prev, "microdecon_water_taxa_prevalence.csv", row.names = FALSE)
cat("\nPrevalence in faeces of the 20 dominant water OTUs (median %):\n")
print(prev %>% group_by(locality) %>%
        summarise(median_prev = median(prev_faeces), max_prev = max(prev_faeces)))

# ── 5. BUILD OBJECTS (host-free, as used in O2/O3) ────────────────────────────
ps_out <- prune_taxa(rownames(otu_c), ps_cur)
otu_table(ps_out) <- otu_table(otu_filt, taxa_are_rows = TRUE)
ps_out <- prune_taxa(taxa_sums(ps_out) > 0, ps_out)
saveRDS(ps_out, paste0("phyloseq_OTU97_LULU_e20_EUK_KKE_curated_", SUFFIX, ".rds"))

ps_out_sp <- tax_glom(ps_out, taxrank = "Species", NArm = FALSE)   # as precuration Step 29
saveRDS(ps_out_sp, paste0("phyloseq_OTU97_LULU_e20_EUK_KKE_curated_", SUFFIX, "_speciesGlom.rds"))
cat("\nSaved", SUFFIX, "| OTUs:", ntaxa(ps_out), "| species-level taxa:", ntaxa(ps_out_sp), "\n")
cat("Faecal samples < 5000 reads after decontamination:\n")
print(sort(colSums(otu_filt[, kk_ids])[colSums(otu_filt[, kk_ids]) < 5000]))

# ── 6. EXPORT FOR SUPPLEMENTARY TABLE S1 ──────────────────────────────────────

w_ids_all <- rownames(sd_w)

s1_otu <- data.frame(
  OTU_ID                     = rownames(otu_c),
  faecal_reads_before        = rowSums(otu_c[, kk_ids, drop = FALSE]),
  faecal_reads_after         = rowSums(otu_filt[, kk_ids, drop = FALSE]),
  faecal_samples_before      = rowSums(otu_c[, kk_ids, drop = FALSE] > 0),
  faecal_samples_after       = rowSums(otu_filt[, kk_ids, drop = FALSE] > 0),
  seawater_reads             = rowSums(otu_w[, w_ids_all, drop = FALSE]),
  seawater_samples           = rowSums(otu_w[, w_ids_all, drop = FALSE] > 0)
) %>%
  mutate(reads_removed_microDecon = faecal_reads_before - faecal_reads_after,
         pct_removed_microDecon   = round(100 * reads_removed_microDecon /
                                            pmax(faecal_reads_before, 1), 1))
write.csv(s1_otu, "S1_microdecon_per_OTU.csv", row.names = FALSE)

# Seawater read counts (all water samples, curated OTUs)
write.csv(data.frame(OTU_ID = rownames(otu_w), otu_w[, w_ids_all, drop = FALSE],
                     check.names = FALSE),
          "S1_seawater_read_counts.csv", row.names = FALSE)

# Faecal read counts after microDecon (the table used in O2/O3)
write.csv(data.frame(OTU_ID = rownames(otu_filt), otu_filt[, kk_ids, drop = FALSE],
                     check.names = FALSE),
          "S1_faecal_read_counts_after_microdecon.csv", row.names = FALSE)

# Water sample metadata (locality and campaign batch)
write.csv(data.frame(sample = w_ids_all, locality = sd_w$locality,
                     batch = w_batch[w_ids_all]),
          "S1_seawater_samples.csv", row.names = FALSE)
cat("\nS1 export written: 4 CSV files (S1_*.csv)\n")
