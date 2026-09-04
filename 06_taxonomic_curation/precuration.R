# ══════════════════════════════════════════════════════════════════════════════
# PRECURATION PIPELINE — METABARCODING DIET SEA URCHINS
# ══════════════════════════════════════════════════════════════════════════════
#
# INPUT:   phyloseq_OTU97_LULU_e20.rds
#
# OUTPUTS:
#   phyloseq_OTU97_LULU_e20_EUK_KKE_curated.rds         ← OTU-level, ALL localities
#   phyloseq_OTU97_LULU_e20_EUK_KKE_curated_speciesGlom.rds  ← species-level ⭐ USE THIS
#   OTUs_removed_log.csv                                 ← all removed OTUs with reason
#   supplementary_OTU_table.csv                          ← full OTU table for SI
#
# PIPELINE STEPS:
#   1  Fix Phaeophyceae missing Phylum
#   2  Keep Eukaryota only
#   3  Change OTUs without Phylum to unknown
#   4  Taxonomic reclassifications (no removal)
#   5  Filter tissue: keep KK + E only (all localities preserved)
#   6  Remove bad samples (AT_2_KK, PC1_58_KK)
#   7  Remove whole phyla: contaminants/non-marine
#   8  Chordata: keep Ascidiacea only
#   9  Arthropoda: remove terrestrial classes
#   10 Amoebozoa: remove all
#   11 Oomycota: remove terrestrial/freshwater genera + Peronosporaceae
#   12 Nematoda: remove terrestrial genera
#   13 Rotifera: remove freshwater genera
#   14 Ciliophora: remove all
#   15 Platyhelminthes: remove Cestoda + Trematoda
#   16 Tardigrada: remove Diploechiniscus (keep Echiniscoides)
#   17 Cnidaria: remove Hydra, Cassiopea, Millepora + tropical Scleractinia;
#               degrade Astrocoeniidae to family level
#   18 Arthropoda: remove freshwater/subterranean genera
#   19 Mollusca: remove freshwater/terrestrial families and orders
#   20 Annelida: remove terrestrial/freshwater Clitellata
#   21 Chlorophyta: remove Chlorella (lab contaminant)
#   22 Remove genera with no marine species
#   23 Remove specific freshwater/terrestrial species
#   24 Remove BLAST-verified bad assignments
#   25 Remove species outside geographic range (EXPANDED — new round)
#   26 Species degradations to genus sp. (EXPANDED — new round)
#   27 Summary + save OTU-level object + removal log
#   28 Build supplementary OTU table (SI)
#   29 Species glom → save final species-level object ⭐
# ══════════════════════════════════════════════════════════════════════════════

library(phyloseq)
library(dplyr)
library(tibble)
library(Biostrings)

setwd(".")  # set to your working directory containing phyloseq_OTU97_LULU_e20.rds

# ── LOAD ──────────────────────────────────────────────────────────────────────
ps <- readRDS("phyloseq_OTU97_LULU_e20.rds")
cat("=== ORIGINAL OBJECT ===\n")
cat("Taxa:", ntaxa(ps), "| Samples:", nsamples(ps), "| Reads:", sum(otu_table(ps)), "\n")

# ── TRACKING: removal log ─────────────────────────────────────────────────────
removed_log <- tibble(
  OTU = character(), reason = character(),
  Domain = character(), Phylum = character(), Class = character(),
  Order = character(), Family = character(), Genus = character(),
  Species = character(), total_reads = numeric()
)

log_removed <- function(ps_before, ps_after, reason_label) {
  gone <- setdiff(taxa_names(ps_before), taxa_names(ps_after))
  if (length(gone) == 0) return(invisible(NULL))
  tax  <- as.data.frame(tax_table(ps_before))[gone, , drop = FALSE]
  otu  <- as(otu_table(ps_before), "matrix")
  if (!taxa_are_rows(ps_before)) otu <- t(otu)
  reads <- rowSums(otu[gone, , drop = FALSE])
  removed_log <<- bind_rows(removed_log, tibble(
    OTU = gone, reason = reason_label,
    Domain = tax$Domain, Phylum = tax$Phylum, Class = tax$Class,
    Order = tax$Order, Family = tax$Family, Genus = tax$Genus,
    Species = tax$Species, total_reads = reads
  ))
  cat("  [-]", length(gone), "OTUs removed —", reason_label, "\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Fix Phaeophyceae with missing Phylum -> Ochrophyta
# ══════════════════════════════════════════════════════════════════════════════
cat("\nSTEP 1 — Fix Phaeophyceae missing Phylum\n")
tax_mat <- as.matrix(tax_table(ps))
tax_mat[!is.na(tax_mat[, "Class"]) &
          tax_mat[, "Class"] == "Phaeophyceae" &
          (is.na(tax_mat[, "Phylum"]) | tax_mat[, "Phylum"] == ""),
        "Phylum"] <- "Ochrophyta"
tax_table(ps) <- tax_table(tax_mat)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Keep Eukaryota only (removes Bacteria ~2.7%, Unassigned ~19.7%)
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 2 — Keep Eukaryota only\n")
ps_before <- ps
ps <- subset_taxa(ps, Domain == "Eukaryota")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Non-Eukaryota (Bacteria / Unassigned domain)")
cat("  Taxa:", ntaxa(ps), "| Reads:", sum(otu_table(ps)), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Remove OTUs without Phylum assignment
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 3 — Remove OTUs without Phylum\n")
ps_before <- ps
ps <- subset_taxa(ps, !is.na(Phylum) & Phylum != "" & Phylum != "Unassigned")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "No Phylum assignment (dinoflagellates, cryptophytes, etc.)")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Taxonomic reclassifications (no OTUs removed)
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 4 — Taxonomic reclassifications\n")
tax_mat <- as.matrix(tax_table(ps))
tax_mat[!is.na(tax_mat[, "Phylum"]) & tax_mat[, "Phylum"] == "Sipuncula",
        "Phylum"] <- "Annelida"       # Sipuncula is subphylum of Annelida
tax_mat[!is.na(tax_mat[, "Phylum"]) & tax_mat[, "Phylum"] == "Discosea",
        "Phylum"] <- "Amoebozoa"      # correct phylum placement
tax_table(ps) <- tax_table(tax_mat)
cat("  Sipuncula -> Annelida | Discosea -> Amoebozoa\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Keep KK (faecal) and E (stomach) tissues only
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 5 — Filter tissue: KK + E (all localities preserved)\n")
ps_before <- ps
ps <- subset_samples(ps, tissue %in% c("KK", "E"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Tissue not KK or E (water/rock samples)")
cat("  Samples:", nsamples(ps), "| Taxa:", ntaxa(ps), "\n")
print(table(sample_data(ps)$tissue, sample_data(ps)$locality))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Remove bad samples
# AT_2_KK: contaminated; PC1_58_KK: sequencing outlier
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 6 — Remove bad samples (AT_2_KK, PC1_58_KK)\n")
ps_before <- ps
ps <- subset_samples(ps, !sample_names(ps) %in% c("AT_2_KK", "PC1_58_KK"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "OTUs lost after removing bad samples (AT_2_KK, PC1_58_KK)")
cat("  Samples:", nsamples(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7 — Remove whole phyla: contaminants and non-marine
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 7 — Remove whole phyla\n")

ps_before <- ps
ps <- subset_taxa(ps, !Phylum %in% c("Imbricatea", "Mucoromycota"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Imbricatea / Mucoromycota: lab contaminants (1-2 reads each)")

ps_before <- ps
ps <- subset_taxa(ps, Phylum != "Hemichordata")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Hemichordata: 8 reads total, water/rock only")

ps_before <- ps
ps <- subset_taxa(ps, !Phylum %in% c("Ascomycota", "Basidiomycota"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Ascomycota / Basidiomycota: terrestrial fungi (Cladosporium, Candida, Mycena, etc.)")

ps_before <- ps
ps <- subset_taxa(ps, Phylum != "Apicomplexa")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Apicomplexa: vertebrate/reptile parasites (Cyclospora, Choleoeimeria, Goussia)")

ps_before <- ps
ps <- subset_taxa(ps, Phylum != "Tubulinea")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Tubulinea: freshwater/terrestrial amoeba (Hartmannella, Hyalosphenia)")

cat("  Taxa after Step 7:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 8 — Chordata: keep Ascidiacea only
# Removed: Actinopteri, Aves, Mammalia, Lepidosauria, Chondrichthyes, NA
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 8 — Chordata: keep Ascidiacea only\n")
ps_before <- ps
ps <- subset_taxa(ps, Phylum != "Chordata" |
                    (!is.na(Class) & Class == "Ascidiacea"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Chordata non-Ascidiacea: fish, birds, mammals, reptiles removed")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 9 — Arthropoda: remove terrestrial classes
# Removed: Insecta, Arachnida, Collembola, Diplopoda, Chilopoda
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 9 — Arthropoda: remove terrestrial classes\n")
ps_before <- ps
ps <- subset_taxa(ps, Phylum != "Arthropoda" |
                    is.na(Class) |
                    !Class %in% c("Insecta", "Arachnida", "Collembola",
                                  "Diplopoda", "Chilopoda"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Arthropoda terrestrial classes: Insecta, Arachnida, Collembola, Diplopoda, Chilopoda")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 10 — Amoebozoa: remove ALL
# Originally Cunea was kept as echinoderm parasite; decision revised.
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 10 — Amoebozoa: remove all\n")
ps_before <- ps
ps <- subset_taxa(ps, Phylum != "Amoebozoa")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Amoebozoa: all removed (Cunea included, decision revised)")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 11 — Oomycota: remove terrestrial/freshwater genera and Peronosporaceae
# Removed: Phytophthora, Hyaloperonospora, Peronospora, Perofascia,
#          Achlya, Saprolegnia, Dictyuchus, Apodachlya, Phytopythium,
#          Pythium (BLAST confirmed terrestrial), Peronosporaceae (family)
# Kept: Aphanomyces (marine parasite), Lagenidium (marine crustacean parasite)
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 11 — Oomycota: remove terrestrial/freshwater\n")
ps_before <- ps
genera_oomy_remove <- c("Phytophthora", "Hyaloperonospora", "Peronospora",
                        "Perofascia", "Achlya", "Saprolegnia",
                        "Dictyuchus", "Apodachlya", "Phytopythium", "Pythium")
ps <- subset_taxa(ps,
                  Phylum != "Oomycota" | (
                    (is.na(Genus)   | !Genus   %in% genera_oomy_remove) &
                      (is.na(Family)  | Family   != "Peronosporaceae")
                  ))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Oomycota terrestrial/freshwater: Phytophthora, Pythium, Saprolegnia, Peronosporaceae, etc.")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 12 — Nematoda: remove terrestrial genera
# Removed: Bursaphelenchus (pine parasite), Panagrellus (soil)
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 12 — Nematoda: remove terrestrial genera\n")
ps_before <- ps
ps <- subset_taxa(ps, Phylum != "Nematoda" |
                    is.na(Genus) |
                    !Genus %in% c("Bursaphelenchus", "Panagrellus"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Nematoda terrestrial genera: Bursaphelenchus (pine), Panagrellus (soil)")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 13 — Rotifera: remove freshwater genera
# Removed: Lecane, Plationus (strictly freshwater),
#          Trichocerca, Brachionus (mostly freshwater)
# Kept: Synchaeta (marine), unclassified
# NEW: Proales fallaciosa removed (freshwater; see Step 25)
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 13 — Rotifera: remove freshwater genera\n")
ps_before <- ps
ps <- subset_taxa(ps, Phylum != "Rotifera" |
                    is.na(Genus) |
                    !Genus %in% c("Lecane", "Plationus", "Trichocerca", "Brachionus"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Rotifera freshwater genera: Lecane, Plationus, Trichocerca, Brachionus")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 14 — Ciliophora: remove all
# Environmental microorganisms; 3 OTUs, 95 reads total
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 14 — Ciliophora: remove all\n")
ps_before <- ps
ps <- subset_taxa(ps, Phylum != "Ciliophora")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Ciliophora: all removed (environmental microorganisms, 95 reads)")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 15 — Platyhelminthes: remove Cestoda and Trematoda
# Vertebrate endoparasites, not sea urchin diet
# Kept: Rhabditophora (free-living flatworms — plausible diet)
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 15 — Platyhelminthes: remove Cestoda + Trematoda\n")
ps_before <- ps
ps <- subset_taxa(ps, Phylum != "Platyhelminthes" |
                    is.na(Class) |
                    !Class %in% c("Cestoda", "Trematoda"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Platyhelminthes Cestoda + Trematoda: vertebrate endoparasites")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 16 — Tardigrada: remove Diploechiniscus (terrestrial/lichen)
# Kept: Echiniscoides (marine intertidal, plausible in diet)
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 16 — Tardigrada: remove Diploechiniscus\n")
ps_before <- ps
ps <- subset_taxa(ps, Phylum != "Tardigrada" |
                    is.na(Genus) |
                    Genus != "Diploechiniscus")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Tardigrada Diploechiniscus: terrestrial/lichen (Echiniscoides kept)")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 17 — Cnidaria: remove freshwater, tropical out-of-range taxa
#
# REMOVE:
#   Hydra (Hydridae): strictly freshwater
#   Cassiopea: Indo-Pacific jellyfish
#   Millepora dichotoma: Indo-Pacific fire coral
#   Scleractinia tropical families: Pocilloporidae, Acroporidae, Poritidae,
#     Siderastreidae; genus Caryophyllia (deep-sea), Montipora (Indo-Pacific)
#   Scleractinia with NA family + NA genus (unclassifiable, very few reads)
#
# DEGRADE to family (no removal):
#   Astrocoeniidae: 1170 reads in 11 samples — spread across samples suggests
#   real signal but species unresolvable -> set Genus + Species to NA
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 17 — Cnidaria: remove freshwater + tropical out-of-range\n")

ps_before <- ps
ps <- subset_taxa(ps, Phylum != "Cnidaria" |
                    is.na(Genus) | Genus != "Hydra")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Cnidaria Hydra: strictly freshwater")

ps_before <- ps
ps <- subset_taxa(ps, Phylum != "Cnidaria" |
                    is.na(Genus) | !Genus %in% c("Cassiopea", "Millepora"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Cnidaria Cassiopea + Millepora: Indo-Pacific out-of-range")

ps_before <- ps
tax_mat <- as.matrix(tax_table(ps))
scleract_remove <- rownames(tax_mat)[
  !is.na(tax_mat[, "Order"]) & tax_mat[, "Order"] == "Scleractinia" & (
    (!is.na(tax_mat[, "Family"]) &
       tax_mat[, "Family"] %in% c("Pocilloporidae", "Acroporidae",
                                  "Poritidae", "Siderastreidae")) |
      (!is.na(tax_mat[, "Genus"]) &
         tax_mat[, "Genus"] %in% c("Caryophyllia", "Montipora")) |
      (is.na(tax_mat[, "Family"]) & is.na(tax_mat[, "Genus"]))
  )
]
ps <- prune_taxa(!taxa_names(ps) %in% scleract_remove, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Scleractinia tropical out-of-range: Pocilloporidae, Acroporidae, Poritidae, Siderastreidae, Caryophyllia, Montipora")

# Degrade Astrocoeniidae to family level (keep OTU, clear genus+species)
tax_mat <- as.matrix(tax_table(ps))
astro <- rownames(tax_mat)[!is.na(tax_mat[, "Family"]) &
                             tax_mat[, "Family"] == "Astrocoeniidae"]
if (length(astro) > 0) {
  tax_mat[astro, "Genus"]   <- NA
  tax_mat[astro, "Species"] <- NA
  tax_table(ps) <- tax_table(tax_mat)
  cat("  [~]", length(astro), "Astrocoeniidae degraded to family level (Genus+Species -> NA)\n")
}
cat("  Taxa after Step 17:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 18 — Arthropoda: remove freshwater/subterranean genera
# Niphargus: subterranean freshwater amphipod
# Alonella: freshwater cladoceran (Chydoridae)
# Moraria: freshwater harpacticoid copepod (Canthocamptidae)
# Amarinus paralacustris: estuarine/freshwater crab
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 18 — Arthropoda: remove freshwater/subterranean genera\n")
ps_before <- ps
tax_mat <- as.matrix(tax_table(ps))
arthro_fw <- rownames(tax_mat)[
  !is.na(tax_mat[, "Phylum"]) & tax_mat[, "Phylum"] == "Arthropoda" & (
    (!is.na(tax_mat[, "Genus"]) &
       tax_mat[, "Genus"] %in% c("Niphargus", "Alonella", "Moraria")) |
      (!is.na(tax_mat[, "Species"]) &
         grepl("Amarinus paralacustris", tax_mat[, "Species"], ignore.case = TRUE))
  )
]
ps <- prune_taxa(!taxa_names(ps) %in% arthro_fw, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Arthropoda freshwater/subterranean: Niphargus, Alonella, Moraria, Amarinus paralacustris")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 19 — Mollusca: remove freshwater and terrestrial taxa
# Families: Planorbidae, Lymnaeidae, Physidae, Thiaridae, Potamididae,
#           Amphibolidae, Assimineidae
# Order: Stylommatophora (all terrestrial snails)
# Genus: Chondrocyclus (terrestrial African gastropod)
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 19 — Mollusca: remove freshwater/terrestrial\n")
ps_before <- ps
tax_mat <- as.matrix(tax_table(ps))
mollusca_remove <- rownames(tax_mat)[
  !is.na(tax_mat[, "Phylum"]) & tax_mat[, "Phylum"] == "Mollusca" & (
    (!is.na(tax_mat[, "Family"]) &
       tax_mat[, "Family"] %in% c("Planorbidae", "Lymnaeidae", "Physidae",
                                  "Thiaridae", "Potamididae", "Amphibolidae",
                                  "Assimineidae")) |
      (!is.na(tax_mat[, "Order"]) & tax_mat[, "Order"] == "Stylommatophora") |
      (!is.na(tax_mat[, "Genus"]) & tax_mat[, "Genus"] == "Chondrocyclus")
  )
]
ps <- prune_taxa(!taxa_names(ps) %in% mollusca_remove, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Mollusca freshwater/terrestrial: Planorbidae, Lymnaeidae, Physidae, Thiaridae, Potamididae, Amphibolidae, Assimineidae, Stylommatophora, Chondrocyclus")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 20 — Annelida: remove terrestrial/freshwater Clitellata
# Crassiclitellata: earthworms (Amynthas, Lumbricidae)
# Limnodriloides, Macquaridrilus: freshwater oligochaetes
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 20 — Annelida: remove terrestrial/freshwater Clitellata\n")
ps_before <- ps
tax_mat <- as.matrix(tax_table(ps))
annelida_remove <- rownames(tax_mat)[
  !is.na(tax_mat[, "Phylum"]) & tax_mat[, "Phylum"] == "Annelida" & (
    (!is.na(tax_mat[, "Order"]) & tax_mat[, "Order"] == "Crassiclitellata") |
      (!is.na(tax_mat[, "Genus"]) &
         tax_mat[, "Genus"] %in% c("Limnodriloides", "Macquaridrilus"))
  )
]
ps <- prune_taxa(!taxa_names(ps) %in% annelida_remove, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Annelida terrestrial/freshwater: Crassiclitellata (Amynthas, Lumbricidae), Limnodriloides, Macquaridrilus")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 21 — Chlorophyta: remove Chlorella (ubiquitous lab contaminant)
# Kept: Micromonas, Bathycoccus (genuine marine picophytoplankton),
#       Codium, Dolichomastix
# NOTE: Bathycoccus + Micromonas also ELIMINATED here per new decision
#       (picoplankton kept in previous pipeline but now removed for diet analyses)
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 21 — Chlorophyta: remove Chlorella + Mamiellophyceae picoplankton\n")
ps_before <- ps
ps <- subset_taxa(ps, Phylum != "Chlorophyta" |
                    is.na(Genus) | Genus != "Chlorella")
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Chlorophyta Chlorella: ubiquitous lab contaminant (Micromonas + Bathycoccus kept)")

# NEW: remove Bathycoccus + Micromonas (Mamiellophyceae — picoplankton,
# not diet; previously kept but now eliminated per curator decision)
ps_before <- ps
ps <- subset_taxa(ps, is.na(Genus) | !Genus %in% c("Bathycoccus", "Micromonas"))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Mamiellophyceae picoplankton: Bathycoccus + Micromonas eliminated (not diet items)")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 22 — Remove genera with no marine species (terrestrial/freshwater)
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 22 — Remove non-marine genera\n")
ps_before <- ps
genera_nonmarine <- c(
  # Annelida freshwater/terrestrial
  "Allolobophora", "Achaeta", "Dero", "Potamothrix", "Dina",
  # Arthropoda terrestrial/freshwater
  "Stenasellus", "Armadillo", "Eucypris",
  # Mollusca terrestrial/freshwater
  "Montenegrina", "Tudorella", "Succinea", "Helix", "Achatinella",
  "Satsuma", "Theba", "Tacheocampylaea", "Trishoplita", "Cochlostoma",
  "Curvella", "Stereophaedusa", "Molema",
  "Bulinus", "Biomphalaria", "Pisidium", "Radix", "Colligyrus",
  "Anisus", "Bellamya", "Pila", "Neoradina", "Parreysia", "Koreosolenaia",
  # Platyhelminthes terrestrial
  "Cephaloflexa", "Luteostriata",
  # Rhodophyta freshwater
  "Kumanoa", "Torularia", "Montagnia", "Virescentia",
  # Bacillariophyta freshwater
  "Gomphonema", "Pinnularia", "Frustulia", "Didymosphenia"
)
tax_df_cur <- as.data.frame(tax_table(ps))
keep_otus <- rownames(tax_df_cur)[is.na(tax_df_cur$Genus) |
                                    !tax_df_cur$Genus %in% genera_nonmarine]
ps <- prune_taxa(keep_otus, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Genera with no marine species (terrestrial/freshwater Annelida, Mollusca, Rhodophyta, Bacillariophyta, etc.)")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 23 — Remove specific freshwater/terrestrial species
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 23 — Remove specific freshwater/terrestrial species\n")
ps_before <- ps
species_nonmarine <- c("Eunotia bilunaris", "Fistulifera pelliculosa",
                       "Stenomelania denisoniensis")
tax_df_cur <- as.data.frame(tax_table(ps))
keep_otus2 <- rownames(tax_df_cur)[is.na(tax_df_cur$Species) |
                                     !tax_df_cur$Species %in% species_nonmarine]
ps <- prune_taxa(keep_otus2, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Specific freshwater/terrestrial species: Eunotia bilunaris, Fistulifera pelliculosa, Stenomelania denisoniensis")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 24 — Remove BLAST-verified bad taxonomic assignments
# Capsalidae: Platyhelminthes, bad assignment
# Oregoniidae: bacterial contamination mislabelled as Arthropoda
# Apolemiidae: diatom mislabelled as Cnidaria
# Pyropia: BLAST confirmed not genuine Pyropia
# Nannochloris: picophytoplankton, environmental contamination
# Mantoniella: picophytoplankton, environmental contamination
# Proboscidactyla: outside geographic range
# Two specific Pyuridae OTUs with no reliable BLAST hit
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 24 — Remove BLAST-verified bad assignments\n")
ps_before <- ps
tax_mat <- as.matrix(tax_table(ps))
blast_remove <- rownames(tax_mat)[
  (!is.na(tax_mat[, "Family"]) & tax_mat[, "Family"] == "Capsalidae")    |
    (!is.na(tax_mat[, "Family"]) & tax_mat[, "Family"] == "Oregoniidae")   |
    (!is.na(tax_mat[, "Family"]) & tax_mat[, "Family"] == "Apolemiidae")   |
    (!is.na(tax_mat[, "Genus"])  & tax_mat[, "Genus"]  == "Pyropia")       |
    (!is.na(tax_mat[, "Genus"])  & tax_mat[, "Genus"]  == "Nannochloris")  |
    (!is.na(tax_mat[, "Genus"])  & tax_mat[, "Genus"]  == "Mantoniella")   |
    (!is.na(tax_mat[, "Genus"])  & tax_mat[, "Genus"]  == "Proboscidactyla") |
    rownames(tax_mat) %in% c("0a55792ca916326e23dfe0a09ad621a5",
                             "99b1600c59520eae5e788c4902f57cee")
]
ps <- prune_taxa(!taxa_names(ps) %in% blast_remove, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "BLAST-verified bad assignments: Capsalidae, Oregoniidae, Apolemiidae, Pyropia, Nannochloris, Mantoniella, Proboscidactyla, Pyuridae no-hit")
cat("  Taxa:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 25 — Remove species outside geographic range (EXPANDED — new round)
#
# NEW additions after second curation audit:
#
# ARTHROPODA:
#   Porcellio spp. — terrestrial isopod (Oniscidea); same group as Spherillo
#   Macrobrachium niphanae — freshwater prawn, SE Asian
#   Attheyella sp. — freshwater harpacticoid copepod
#   Gammarus minus / G. lacustris — freshwater amphipods; KEEP 2 OTUs at
#     family level (Gammaridae): a98da7992e83dd20bbfcac8178d16534,
#     bd0bf6019297e578ecc7fc44d8460a0e; rest of Gammarus OTUs eliminated
#   Sinorchestia taiwanensis — semiterrestrial Asian amphipod (degraded Step 26)
#   Callinectes spp. — W Atlantic portunid; BLAST confirmed: ELIMINATE ALL
#   Lophopanopeus bellus — NE Pacific crab, BLAST artefact (degraded Step 26)
#   Pugettia spp. — NE Pacific kelp crabs (degraded Step 26)
#   Xeruca formosensis — Indo-Pacific Ocypodidae (degraded Step 26)
#
# OCHROPHYTA:
#   Hormosira banksii — endemic Australia/NZ; BLAST confirmed: degrade to Fucales
#   Pelagophycus porra — California kelp; BLAST confirmed: degrade to Fucales
#   Durvillaea spp. — Southern Ocean; BLAST confirmed: degrade to Fucales
#   Pterygophora californica — NE Pacific; direct removal
#   Egregia menziesii — NE Pacific; direct removal
#   Costaria costata — N Pacific; direct removal
#
# RHODOPHYTA (BLAST verified — ALL KEPT, degraded to Order level):
#   Mychodeophyllum papillitectum — Southern Ocean/Australia; direct removal
#   Mychodea spp. — Australian genus; BLAST: keep degraded to Order
#   Bossiella spp. — NE Pacific corallines; BLAST: keep degraded to Order
#   Acanthopeltis japonica — Japan/NW Pacific; BLAST: keep degraded to Family
#   Gloiopeltis tenax — East Asia; BLAST: keep degraded to Order
#   Austrokallymenia rosea — Southern Ocean; direct removal
#   Entwisleia bella — Australian endemic; direct removal
#   Fuscifolium papenfussii — Pacific; direct removal
#   Titanophycus validus / Titanophora weberae — Indo-Pacific; direct removal
#   Stauromenia sp. — Indo-Pacific/Australian; BLAST: keep degraded to Order
#   Meristotheca sp. — Indo-Pacific; BLAST: keep degraded to Order
#   Faucheocolax attenuata — Pacific parasite; BLAST: keep degraded to Order
#
# ANNELIDA:
#   Tubifex tubifex — freshwater oligochaete (degraded Step 26)
#   Pterobdella sp. — fish leech; direct removal
#   Clitellata — BLAST verified: REMOVE ALL except OTU
#     86102227f65ad9050ca97a771a11095e (= Lumbrineridae, Polychaeta)
#     which is reclassified to Polychaeta in Step 26
#
# MOLLUSCA:
#   Lunella spp. — Indo-Pacific Turbinidae (degraded Step 26)
#   Uroteuthis sibogae — Indo-Pacific squid; direct removal
#   Indothais rufotincta — Indo-Pacific Muricidae (degraded Step 26)
#   Californiconus californicus — NE Pacific; direct removal
#   Ocinebrellus inornatus — NE Pacific muricid; direct removal
#   Omphalius nigerrimus — NW Pacific trochid; direct removal
#   Neritina canalis — Indo-Pacific freshwater nerite (degraded Step 26)
#   Geukensia granosissima — W Atlantic mussel; direct removal
#   Graneledone taniwha — deep-sea NZ octopus; BLAST confirmed: ELIMINATE
#   Lepetodrilus pustulosus — hydrothermal vent limpet (degraded Step 26)
#   Strubellia sp. — Indo-Pacific acochlidian (degraded Step 26)
#   Naytia glabrata — Indo-Pacific (degraded Step 26)
#   Montfortula variegata — Indo-Pacific fissurellid (degraded Step 26)
#
# CNIDARIA:
#   Catostylus sp. — Indo-Pacific jellyfish; BLAST confirmed: ELIMINATE ALL
#   Phyllorhiza punctata — Indo-Pacific jellyfish; direct removal
#   Paraminabea aldersladei — Indo-Pacific soft coral; BLAST confirmed: ELIMINATE
#
# ECHINODERMATA:
#   Linckia laevigata — Indo-Pacific starfish; direct removal
#   Cenolia spp. — Indo-Pacific crinoid; direct removal
#   Comanthus sp. — Indo-Pacific crinoid; direct removal
#   Florometra sp. — N Pacific crinoid; direct removal
#   Mellita sp. — W Atlantic sand dollar; direct removal
#   Arbacia spatuligera — SE Pacific urchin; direct removal
#
# NEMERTEA:
#   Amphiporus texanus — NW Atlantic/Gulf (degraded to genus in Step 26)
#
# ROTIFERA:
#   Proales fallaciosa — freshwater; direct removal
#
# CTENOPHORA:
#   Coeloplana loyai — benthic Indo-Pacific ctenophore; direct removal
#
# TAXA CONSERVED (no action):
#   Erenna laciniata — cosmopolitan mesopelagic siphonophore (degraded Step 26)
#   Eurythenes gryllus — cosmopolitan abyssal amphipod (degraded Step 26)
#   Percnon gibbesi — invasive crab established in Canaries since 1999 (KEPT)
#   Benthodytes sanguinolenta — Atlantic holothurian (degraded Step 26)
#   Laminaria / Saccharina spp. — possible diet items (degraded Step 26)
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 25 — Remove out-of-range species (expanded second round)\n")

# ---- 25a: DIRECT REMOVALS (genera fully eliminated) -------------------------
ps_before <- ps
tax_mat <- as.matrix(tax_table(ps))

genera_oor_remove <- c(
  # Arthropoda
  "Porcellio",          # terrestrial isopod
  "Macrobrachium",      # freshwater prawn SE Asia
  "Attheyella",         # freshwater harpacticoid
  "Callinectes",        # W Atlantic portunid — BLAST confirmed: no Canary hits
  # Ochrophyta
  "Pterygophora",       # NE Pacific kelp
  "Egregia",            # NE Pacific
  "Costaria",           # N Pacific
  # Rhodophyta
  "Mychodeophyllum",    # Southern Ocean/Australia
  "Austrokallymenia",   # Southern Ocean
  "Entwisleia",         # Australia endemic
  "Fuscifolium",        # Pacific
  "Titanophycus",       # Indo-Pacific
  "Titanophora",        # Indo-Pacific
  # Annelida
  "Pterobdella",        # fish leech (parasite, not diet)
  # Mollusca
  "Uroteuthis",         # Indo-Pacific squid
  "Californiconus",     # NE Pacific
  "Ocinebrellus",       # NE Pacific muricid
  "Omphalius",          # NW Pacific trochid
  "Geukensia",          # W Atlantic mussel
  "Graneledone",        # deep-sea NZ octopus — BLAST confirmed: eliminate
  # Cnidaria
  "Phyllorhiza",        # Indo-Pacific jellyfish
  "Catostylus",         # Indo-Pacific jellyfish — BLAST confirmed: eliminate
  "Paraminabea",        # Indo-Pacific soft coral — BLAST confirmed: eliminate
  # Echinodermata
  "Linckia",            # Indo-Pacific starfish
  "Cenolia",            # Indo-Pacific crinoid
  "Comanthus",          # Indo-Pacific crinoid
  "Florometra",         # N Pacific crinoid
  "Mellita"             # W Atlantic sand dollar
)

species_oor_remove <- c(
  "Arbacia spatuligera"  # SE Pacific urchin — mislabelled A. lixula
)

oor_new <- rownames(tax_mat)[
  (!is.na(tax_mat[, "Genus"])   & tax_mat[, "Genus"]   %in% genera_oor_remove) |
    (!is.na(tax_mat[, "Species"]) & tax_mat[, "Species"] %in% species_oor_remove)
]
ps <- prune_taxa(!taxa_names(ps) %in% oor_new, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps,
            "Out-of-range genera/species (new round, BLAST verified): Porcellio, Macrobrachium, Attheyella, Callinectes, Pterygophora, Egregia, Costaria, Mychodeophyllum, Austrokallymenia, Entwisleia, Fuscifolium, Titanophycus, Titanophora, Pterobdella, Uroteuthis, Californiconus, Ocinebrellus, Omphalius, Geukensia, Graneledone, Phyllorhiza, Catostylus, Paraminabea, Linckia, Cenolia, Comanthus, Florometra, Mellita, Arbacia spatuligera")

# ---- 25b: DIRECT REMOVALS (specific species + targeted OTU-level) -----------
ps_before <- ps
species_oor_2 <- c(
  "Proales fallaciosa",   # freshwater rotifer
  "Coeloplana loyai"      # Indo-Pacific benthic ctenophore
)
tax_df_c <- as.data.frame(tax_table(ps))
keep_oor2 <- rownames(tax_df_c)[is.na(tax_df_c$Species) |
                                  !tax_df_c$Species %in% species_oor_2]
ps <- prune_taxa(keep_oor2, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps,
            "Out-of-range species (new round 2): Proales fallaciosa (freshwater rotifer), Coeloplana loyai (Indo-Pacific ctenophore)")

# ---- 25c: Clitellata — remove ALL except one BLAST-verified Polychaeta OTU --
# OTU 86102227f65ad9050ca97a771a11095e is Lumbrineridae (Polychaeta), not
# Clitellata — it was misclassified. It is KEPT and reclassified in Step 26.
# All other Clitellata OTUs (Tubifex, Pterobdella already gone; unassigned ones)
# are removed here.
ps_before <- ps
tax_df_c2 <- as.data.frame(tax_table(ps))
lumbrinerid_otu <- "86102227f65ad9050ca97a771a11095e"
clitellata_remove <- rownames(tax_df_c2)[
  !is.na(tax_df_c2$Class) &
    tax_df_c2$Class == "Clitellata" &
    rownames(tax_df_c2) != lumbrinerid_otu
]
ps <- prune_taxa(!taxa_names(ps) %in% clitellata_remove, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps,
            "Clitellata: all remaining OTUs removed (freshwater/terrestrial); Lumbrineridae OTU 86102227f65ad9050ca97a771a11095e kept and reclassified to Polychaeta")

# ---- 25d: Gammarus — remove all EXCEPT two BLAST-verified family-level OTUs --
# OTUs a98da7992e83dd20bbfcac8178d16534 and bd0bf6019297e578ecc7fc44d8460a0e
# BLAST to marine/coastal Gammaridae; kept at family level (degraded in Step 26).
# All other Gammarus OTUs are eliminated.
ps_before <- ps
tax_df_c3 <- as.data.frame(tax_table(ps))
gammarus_keep <- c("a98da7992e83dd20bbfcac8178d16534",
                   "bd0bf6019297e578ecc7fc44d8460a0e")
gammarus_remove <- rownames(tax_df_c3)[
  !is.na(tax_df_c3$Genus) &
    tax_df_c3$Genus == "Gammarus" &
    !rownames(tax_df_c3) %in% gammarus_keep
]
ps <- prune_taxa(!taxa_names(ps) %in% gammarus_remove, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps,
            "Gammarus: OTUs with no Atlantic hit removed; 2 OTUs kept at Gammaridae family level (BLAST verified)")

# ---- 25e: original out-of-range block (kept from Step 25 first round) -------
ps_before <- ps
tax_mat <- as.matrix(tax_table(ps))
oor_genera_orig  <- c("Sulcospira", "Stenothyra", "Wallaconchis", "Platevindex")
oor_species_orig <- c("Sepiadarium kochi", "Euprymna berryi", "Alderia willowi")
oor_remove_orig  <- rownames(tax_mat)[
  (!is.na(tax_mat[, "Genus"])   & tax_mat[, "Genus"]   %in% oor_genera_orig)  |
    (!is.na(tax_mat[, "Species"]) & tax_mat[, "Species"] %in% oor_species_orig)
]
ps <- prune_taxa(!taxa_names(ps) %in% oor_remove_orig, ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
log_removed(ps_before, ps, "Out-of-range species: Indo-Pacific / freshwater Asian / Eastern Pacific")
cat("  Taxa after Step 25:", ntaxa(ps), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 26 — Species degradations to genus/family/class sp.
#           (OTUs KEPT; only Species and sometimes Genus field cleared)
#
# NEW additions (second curation round):
#   - Out-of-range genera/species degraded to higher rank rather than
#     fully removed, to preserve the read signal at coarser resolution.
#   - "DEGRADE TO ORDER": Species + Genus cleared (assigned at Order level)
#   - "DEGRADE TO FAMILY": Species + Genus cleared (assigned at Family level)
#   - "DEGRADE TO CLASS": Species + Genus + Family + Order cleared
#   - "DEGRADE TO GENUS sp.": Species field replaced with "Genus sp."
# ══════════════════════════════════════════════════════════════════════════════
cat("STEP 26 — Species degradations (genus sp.) — expanded\n")
tax_mat <- as.matrix(tax_table(ps))

# Helper: degrade Species only (genus sp.)
degrade_species <- function(tax_mat, sp_pattern, new_sp) {
  idx <- !is.na(tax_mat[, "Species"]) &
    grepl(sp_pattern, tax_mat[, "Species"], ignore.case = TRUE)
  if (any(idx)) tax_mat[idx, "Species"] <- new_sp
  tax_mat
}

# Helper: degrade to Order level (clear Genus + Species)
degrade_to_order <- function(tax_mat, genus_names) {
  idx <- !is.na(tax_mat[, "Genus"]) & tax_mat[, "Genus"] %in% genus_names
  if (any(idx)) {
    tax_mat[idx, "Genus"]   <- NA
    tax_mat[idx, "Species"] <- NA
  }
  tax_mat
}

# Helper: degrade to Family level (clear Genus + Species)
degrade_to_family <- function(tax_mat, genus_names) {
  idx <- !is.na(tax_mat[, "Genus"]) & tax_mat[, "Genus"] %in% genus_names
  if (any(idx)) {
    tax_mat[idx, "Genus"]   <- NA
    tax_mat[idx, "Species"] <- NA
  }
  tax_mat
}

# Helper: degrade to Class level (clear Order + Family + Genus + Species)
degrade_to_class <- function(tax_mat, genus_names) {
  idx <- !is.na(tax_mat[, "Genus"]) & tax_mat[, "Genus"] %in% genus_names
  if (any(idx)) {
    tax_mat[idx, "Order"]   <- NA
    tax_mat[idx, "Family"]  <- NA
    tax_mat[idx, "Genus"]   <- NA
    tax_mat[idx, "Species"] <- NA
  }
  tax_mat
}

# ── ARTHROPODA ────────────────────────────────────────────────────────────────
# Callinectes: eliminated in Step 25 — no degradation needed here.

# DEGRADE TO FAMILY — Gammaridae: only 2 surviving OTUs (rest removed Step 25d)
# The two kept OTUs already have Genus=Gammarus; clear to family level.
gammarus_keep_ids <- c("a98da7992e83dd20bbfcac8178d16534",
                       "bd0bf6019297e578ecc7fc44d8460a0e")
idx_gam <- rownames(tax_mat) %in% gammarus_keep_ids
if (any(idx_gam)) {
  tax_mat[idx_gam, "Genus"]   <- NA
  tax_mat[idx_gam, "Species"] <- NA
}

# DEGRADE TO CLASS — Malacostraca (Sinorchestia: semiterrestrial Asian genus)
tax_mat <- degrade_to_class(tax_mat, c("Sinorchestia"))

# DEGRADE TO FAMILY — Epialtidae (Pugettia: NE Pacific kelp crabs)
tax_mat <- degrade_to_family(tax_mat, c("Pugettia", "Lophopanopeus"))

# DEGRADE TO GENUS sp. — Xeruca (Ocypodidae Indo-Pacific: degrade species only)
tax_mat <- degrade_species(tax_mat, "Xeruca formosensis", "Xeruca sp.")

# ── OCHROPHYTA ────────────────────────────────────────────────────────────────
# BLAST result: Hormosira, Pelagophycus and Durvillaea best hits are Fucales
# (Sargassaceae / Fucaceae relatives present in Canaries/Macaronesia).
# Degrade to ORDER = Fucales (clear Genus + Species only).
tax_mat <- degrade_to_order(tax_mat, c(
  "Hormosira",    # endemic AUS/NZ -> BLAST: Fucales
  "Pelagophycus", # California kelp -> BLAST: Fucales
  "Durvillaea"    # Southern Ocean -> BLAST: Fucales
))

# ── RHODOPHYTA ────────────────────────────────────────────────────────────────
# DEGRADE TO ORDER — clear Genus + Species (keep Order-level signal)
tax_mat <- degrade_to_order(tax_mat, c(
  "Mychodea",       # Australian genus -> Gigartinales
  "Bossiella",      # NE Pacific coralline -> Corallinales
  "Gloiopeltis",    # East Asia -> Gigartinales
  "Stauromenia",    # Indo-Pacific -> Gigartinales
  "Meristotheca",   # Indo-Pacific -> Gigartinales
  "Faucheocolax"    # Pacific parasite rhodophyte -> Rhodymeniales
))

# DEGRADE TO FAMILY — clear Genus + Species
tax_mat <- degrade_to_family(tax_mat, c(
  "Acanthopeltis"   # Japan/NW Pacific Gelidiales -> Gelidiaceae
))

# ── ANNELIDA ─────────────────────────────────────────────────────────────────
# RECLASSIFY — OTU 86102227f65ad9050ca97a771a11095e was assigned to Clitellata
# but BLAST confirms it is Lumbrineridae (Polychaeta: Eunicida). Fix taxonomy.
lumbrinerid_id <- "86102227f65ad9050ca97a771a11095e"
if (lumbrinerid_id %in% rownames(tax_mat)) {
  tax_mat[lumbrinerid_id, "Class"]  <- "Polychaeta"
  tax_mat[lumbrinerid_id, "Order"]  <- "Eunicida"
  tax_mat[lumbrinerid_id, "Family"] <- "Lumbrineridae"
  tax_mat[lumbrinerid_id, "Genus"]  <- NA
  tax_mat[lumbrinerid_id, "Species"] <- NA
  cat("  [~] OTU 86102227f65ad9050ca97a771a11095e reclassified: Clitellata -> Polychaeta / Lumbrineridae\n")
}

# DEGRADE TO FAMILY — Tubifex (freshwater oligochaete; keep Naididae signal)
# Note: most Clitellata already eliminated in Step 25c; this is a safeguard.
tax_mat <- degrade_to_family(tax_mat, c("Tubifex"))

# DEGRADE TO CLASS (Clitellata) — any unassigned Clitellata OTUs that survived.
# No action needed — they will remain as Clitellata sp. after glom.

# ── MOLLUSCA ─────────────────────────────────────────────────────────────────
# DEGRADE TO FAMILY — Turbinidae (Lunella: Indo-Pacific)
tax_mat <- degrade_to_family(tax_mat, c("Lunella"))

# DEGRADE TO FAMILY — Muricidae (Indothais: Indo-Pacific)
tax_mat <- degrade_to_family(tax_mat, c("Indothais"))

# DEGRADE TO GASTROPODA (class) — Neritina (freshwater Indo-Pacific nerite)
tax_mat <- degrade_to_class(tax_mat, c("Neritina"))

# Graneledone: eliminated in Step 25 — no degradation needed.

# DEGRADE TO FAMILY — Lepetodrilidae (vent limpet; keep family signal)
tax_mat <- degrade_to_family(tax_mat, c("Lepetodrilus"))

# DEGRADE TO FAMILY — for remaining Indo-Pacific molluscs
tax_mat <- degrade_to_family(tax_mat, c(
  "Strubellia",    # acochlidian Indo-Pacific
  "Naytia",        # Indo-Pacific nassariid
  "Montfortula"    # Indo-Pacific fissurellid
))

# ── CNIDARIA ─────────────────────────────────────────────────────────────────
# Catostylus: eliminated in Step 25 — no degradation needed.
# Paraminabea: eliminated in Step 25 — no degradation needed.

# ── ECHINODERMATA ─────────────────────────────────────────────────────────────
# None remaining — all out-of-range echinoderms removed in Step 25.

# ── NEMERTEA ──────────────────────────────────────────────────────────────────
# DEGRADE TO GENUS sp. — Amphiporus texanus (170k reads; NW Atlantic species
# but Amphiporus genus present in Atlantic; keep genus-level signal)
tax_mat <- degrade_species(tax_mat, "Amphiporus texanus", "Amphiporus sp.")

# ── TAXA CONSERVED WITH DEGRADATION ──────────────────────────────────────────
# Erenna laciniata -> degrade to Class Hydrozoa (cosmopolitan deep siphonophore)
tax_mat <- degrade_to_class(tax_mat, c("Erenna"))

# Eurythenes gryllus -> degrade to Family Lysianassidae (cosmopolitan abyssal amphipod)
tax_mat <- degrade_to_family(tax_mat, c("Eurythenes"))

# Benthodytes sanguinolenta -> degrade to Family Psychropotidae (Atlantic holothurian)
tax_mat <- degrade_to_family(tax_mat, c("Benthodytes"))

# Laminaria / Saccharina -> degrade to Family Laminariaceae (possible diet items)
tax_mat <- degrade_to_family(tax_mat, c("Laminaria", "Saccharina"))

# ── ORIGINAL DEGRADATIONS (from Step 26 first round) ─────────────────────────
degrade_map_orig <- list(
  # Nemertea
  list(sp = "Berkeleya fennica",        new = "Berkeleya sp."),
  list(sp = "Hemiaulus sinensis",       new = "Hemiaulus sp."),
  list(sp = "Navicula glaciei",         new = "Navicula sp."),
  list(sp = "Navicula cryptocephala",   new = "Navicula sp."),
  list(sp = "Sellaphora pupula",        new = "Sellaphora sp."),
  list(sp = "Nitzschia palea",          new = "Nitzschia sp."),
  list(sp = "Lomentaria orcadensis",    new = "Lomentaria sp."),
  list(sp = "Hormosira banksii",        new = "Hormosira sp."),
  list(sp = "Coccophora langsdorfii",   new = "Coccophora sp."),
  list(sp = "Pelagophycus porra",       new = "Pelagophycus sp."),
  list(sp = "Laminaria solidungula",    new = "Laminaria sp."),
  list(sp = "Phasianella australis",    new = "Phasianella sp."),
  list(sp = "Haliotis rubra",           new = "Haliotis sp."),
  list(sp = "Nerita atramentosa",       new = "Nerita sp."),
  list(sp = "Neritina canalis",         new = "Neritina sp."),
  list(sp = "Merica sinensis",          new = "Merica sp."),
  list(sp = "Gammarus lacustris",       new = "Gammarus sp."),
  list(sp = "Gammarus minus",           new = "Gammarus sp."),
  list(sp = "Alpheus bellulus",         new = "Alpheus sp."),
  list(sp = "Alpheus dolerus",          new = "Alpheus sp."),
  list(sp = "Penaeus merguiensis",      new = "Penaeus sp."),
  list(sp = "Pugettia producta",        new = "Pugettia sp."),
  list(sp = "Lophopanopeus bellus",     new = "Lophopanopeus sp."),
  list(sp = "Sinorchestia taiwanensis", new = "Sinorchestia sp."),
  list(sp = "Xeruca formosensis",       new = "Xeruca sp."),
  list(sp = "Ampithoe dalli",           new = "Ampithoe sp."),
  list(sp = "Macrobrachium niphanae",   new = "Macrobrachium sp."),
  list(sp = "Oratosquilla oratoria",    new = "Oratosquilla sp."),
  list(sp = "Sicyonia brevirostris",    new = "Sicyonia sp."),
  list(sp = "Solenocera crassicornis",  new = "Solenocera sp."),
  list(sp = "Xanthias latifrons",       new = "Xanthias sp."),
  list(sp = "Eucrate crenata",          new = "Eucrate sp."),
  list(sp = "Charybdis erythrodactyla", new = "Charybdis sp."),
  list(sp = "Puerulus sewelli",         new = "Puerulus sp."),
  list(sp = "Cloridopsis scorpio",      new = "Cloridopsis sp."),
  list(sp = "Porcellio spinicornis",    new = "Porcellio sp."),
  list(sp = "Hydroides ezoensis",          new = "Hydroides sp."),
  list(sp = "Apionsoma misakianum",        new = "Apionsoma sp."),
  list(sp = "Octobranchus lingulatus",     new = "Octobranchus sp."),
  list(sp = "Tubifex tubifex",             new = "Tubifex sp."),
  list(sp = "Ancistrosyllis groenlandica", new = "Ancistrosyllis sp."),
  list(sp = "Anadara broughtonii",      new = "Anadara sp."),
  list(sp = "Conus ebraeus",            new = "Conus sp."),
  list(sp = "Wallaconchis buetschlii",  new = "Wallaconchis sp."),
  list(sp = "Platevindex luteus",       new = "Platevindex sp."),
  list(sp = "Platevindex applanatus",   new = "Platevindex sp."),
  list(sp = "Callyspongia vaginalis",   new = "Callyspongia sp."),
  list(sp = "Securiflustra securifrons", new = "Securiflustra sp."),
  list(sp = "Nanomia cara",             new = "Nanomia sp.")
)

n_degraded <- 0
for (item in degrade_map_orig) {
  idx <- !is.na(tax_mat[, "Species"]) &
    grepl(item$sp, tax_mat[, "Species"], ignore.case = TRUE)
  if (any(idx)) {
    tax_mat[idx, "Species"] <- item$new
    n_degraded <- n_degraded + sum(idx)
  }
}
tax_table(ps) <- tax_table(tax_mat)
cat("  [~] Species degradations applied (genus sp. + order/family/class level)\n")
cat("  [~]", n_degraded, "additional OTUs degraded via species pattern matching\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 27 — Summary + save OTU-level curated object + removal log
# ══════════════════════════════════════════════════════════════════════════════
cat("\n╔══════════════════════════════════════════════════════╗\n")
cat("║   CURATED OBJECT — OTU LEVEL (all localities)       ║\n")
cat("╚══════════════════════════════════════════════════════╝\n")
cat("Samples:", nsamples(ps), "\n")
cat("Taxa (OTUs):", ntaxa(ps), "\n")
cat("Total reads:", sum(otu_table(ps)), "\n")
cat("\nSamples by tissue and locality:\n")
print(table(sample_data(ps)$tissue, sample_data(ps)$locality))

tax_df_s <- as.data.frame(tax_table(ps))
otu_mat_s <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu_mat_s <- t(otu_mat_s)
cat("\nPhylum composition:\n")
data.frame(Phylum = tax_df_s$Phylum, reads = rowSums(otu_mat_s)) %>%
  group_by(Phylum) %>%
  summarise(n_OTUs = n(), total_reads = sum(reads), .groups = "drop") %>%
  mutate(pct_reads = round(100 * total_reads / sum(total_reads), 2)) %>%
  arrange(desc(total_reads)) %>%
  print(n = 40)

saveRDS(ps, "phyloseq_OTU97_LULU_e20_EUK_KKE_curated.rds")
cat("\nSaved: phyloseq_OTU97_LULU_e20_EUK_KKE_curated.rds\n")

write.csv(removed_log, "OTUs_removed_log.csv", row.names = FALSE)
cat("Saved: OTUs_removed_log.csv (", nrow(removed_log), "OTUs logged)\n")

cat("\nRemoval summary by reason:\n")
removed_log %>%
  group_by(reason) %>%
  summarise(n_OTUs = n(), total_reads = sum(total_reads), .groups = "drop") %>%
  arrange(desc(n_OTUs)) %>%
  print(n = 40)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 28 — Build supplementary OTU table for SI
# Contains ALL OTUs (retained + removed), with columns:
#   OTU_ID | Domain-Species taxonomy | total_reads | pct_reads |
#   prevalence_nSamps | removed (YES/no) | removal_reason |
#   sequence | reads per sample (wide)
# OTUs marked removed = YES should be greyed out in the SI table
# ══════════════════════════════════════════════════════════════════════════════
cat("\nSTEP 28 — Building supplementary OTU table...\n")

# Rebuild minimal pre-removal base (Steps 1-6 only) as reference universe
ps_si <- readRDS("phyloseq_OTU97_LULU_e20.rds")
tax_si_mat <- as.matrix(tax_table(ps_si))
tax_si_mat[!is.na(tax_si_mat[, "Class"]) &
             tax_si_mat[, "Class"] == "Phaeophyceae" &
             (is.na(tax_si_mat[, "Phylum"]) | tax_si_mat[, "Phylum"] == ""),
           "Phylum"] <- "Ochrophyta"
tax_table(ps_si) <- tax_table(tax_si_mat)
ps_si <- subset_taxa(ps_si, Domain == "Eukaryota")
ps_si <- subset_taxa(ps_si, !is.na(Phylum) & Phylum != "" & Phylum != "Unassigned")
tax_si_mat <- as.matrix(tax_table(ps_si))
tax_si_mat[!is.na(tax_si_mat[, "Phylum"]) &
             tax_si_mat[, "Phylum"] == "Sipuncula", "Phylum"] <- "Annelida"
tax_si_mat[!is.na(tax_si_mat[, "Phylum"]) &
             tax_si_mat[, "Phylum"] == "Discosea",  "Phylum"] <- "Amoebozoa"
tax_table(ps_si) <- tax_table(tax_si_mat)
ps_si <- subset_samples(ps_si, tissue %in% c("KK", "E"))
ps_si <- subset_samples(ps_si, !sample_names(ps_si) %in% c("AT_2_KK", "PC1_58_KK"))
ps_si <- prune_taxa(taxa_sums(ps_si) > 0, ps_si)

tax_si  <- as.data.frame(tax_table(ps_si))
otu_si  <- as(otu_table(ps_si), "matrix")
if (!taxa_are_rows(ps_si)) otu_si <- t(otu_si)

has_ref <- !is.null(tryCatch(refseq(ps_si), error = function(e) NULL))
seq_vec <- if (has_ref) as.character(refseq(ps_si))[rownames(tax_si)] else
  rep(NA_character_, nrow(tax_si))

si_table <- tax_si %>%
  rownames_to_column("OTU_ID") %>%
  mutate(
    total_reads       = rowSums(otu_si[OTU_ID, , drop = FALSE]),
    pct_reads         = round(100 * total_reads / sum(otu_si), 5),
    prevalence_nSamps = rowSums(otu_si[OTU_ID, , drop = FALSE] > 0),
    removed           = ifelse(OTU_ID %in% removed_log$OTU, "YES", "no"),
    removal_reason    = removed_log$reason[match(OTU_ID, removed_log$OTU)],
    sequence          = seq_vec[OTU_ID]
  ) %>%
  arrange(Phylum, desc(total_reads))

sample_read_cols <- as.data.frame(otu_si[si_table$OTU_ID, , drop = FALSE])
si_table <- bind_cols(si_table, sample_read_cols)

write.csv(si_table, "supplementary_OTU_table.csv", row.names = FALSE)
cat("Saved: supplementary_OTU_table.csv\n")
cat("  Total OTUs:", nrow(si_table),
    "| Removed:", sum(si_table$removed == "YES"),
    "| Retained:", sum(si_table$removed == "no"), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 29 — Species glom → FINAL OBJECT FOR ANALYSES
# tax_glom collapses OTUs sharing the same species assignment.
# NArm = FALSE preserves OTUs without species (assigned to higher ranks).
# ══════════════════════════════════════════════════════════════════════════════
cat("\nSTEP 29 — Species glom\n")
cat("  Before glom — Taxa (OTUs):", ntaxa(ps), "| Reads:", sum(otu_table(ps)), "\n")

ps_final <- tax_glom(ps, taxrank = "Species", NArm = FALSE)

tax_sp <- as.data.frame(tax_table(ps_final))
cat("  After glom  — Taxa:", ntaxa(ps_final), "| Reads:", sum(otu_table(ps_final)), "\n")
cat("  Species assigned:", sum(!is.na(tax_sp$Species) & tax_sp$Species != ""), "\n")
cat("  Species NA:      ", sum(is.na(tax_sp$Species) | tax_sp$Species == ""), "\n")

saveRDS(ps_final, "phyloseq_OTU97_LULU_e20_EUK_KKE_curated_speciesGlom.rds")

# ══════════════════════════════════════════════════════════════════════════════
# DONE
# ══════════════════════════════════════════════════════════════════════════════
cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║   PRECURATION COMPLETE                                       ║\n")
cat("╠══════════════════════════════════════════════════════════════╣\n")
cat("║ DATASETS FOR ANALYSES:                                       ║\n")
cat("║                                                              ║\n")
cat("║  ⭐ OBJECTIVE 1 (KK vs E, Blanes):                          ║\n")
cat("║     phyloseq_OTU97_LULU_e20_EUK_KKE_curated_speciesGlom.rds ║\n")
cat("║     -> subset_samples(ps, locality == 'Blanes')              ║\n")
cat("║                                                              ║\n")
cat("║  ⭐ OBJECTIVE 2 (Vent vs Control, La Palma):                 ║\n")
cat("║     phyloseq_OTU97_LULU_e20_EUK_KKE_curated_speciesGlom.rds ║\n")
cat("║     -> subset_samples(ps, locality != 'Blanes')              ║\n")
cat("║                                                              ║\n")
cat("║  OTU-level (if needed):                                      ║\n")
cat("║     phyloseq_OTU97_LULU_e20_EUK_KKE_curated.rds             ║\n")
cat("║                                                              ║\n")
cat("║  Supplementary table SI:  supplementary_OTU_table.csv        ║\n")
cat("║  Removed OTUs log:        OTUs_removed_log.csv               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")


# ══════════════════════════════════════════════════════════════════════════════
# BONUS: Extract sequences for manual BLAST verification
#        (run this block AFTER the main pipeline, on the CURATED object)
#
# All taxa below have been BLAST-verified and decisions applied above.
# This block is kept as a record and for any future re-verification.
#
# ══════════════════════════════════════════════════════════════════════════════

# ── LOAD the curated OTU-level object (pre-glom) ─────────────────────────────
ps_cur <- readRDS("phyloseq_OTU97_LULU_e20_EUK_KKE_curated.rds")

# ── Extract reference sequences ───────────────────────────────────────────────
ref_all <- refseq(ps_cur)
tax_cur <- as.data.frame(tax_table(ps_cur))

# Helper: extract FASTA for a vector of OTU IDs
write_fasta <- function(otu_ids, filename, tax_df, ref_seqs) {
  ids_present <- intersect(otu_ids, names(ref_seqs))
  if (length(ids_present) == 0) {
    cat("  [!] No sequences found for:", filename, "\n")
    return(invisible(NULL))
  }
  seqs <- ref_seqs[ids_present]
  # Add taxonomy to FASTA header for readability
  headers <- sapply(ids_present, function(id) {
    tx <- tax_df[id, ]
    paste0(id, " | ",
           paste(na.omit(c(tx$Phylum, tx$Class, tx$Order,
                           tx$Family, tx$Genus, tx$Species)),
                 collapse = "; "))
  })
  names(seqs) <- headers
  writeXStringSet(seqs, filepath = filename, format = "fasta")
  cat("  [+] Written:", filename, "(", length(seqs), "sequences )\n")
}

dir.create("sequences_to_verify", showWarnings = FALSE)

# ── 1. Mollusca — still pending verification ──────────────────────────────────
# (Graneledone already eliminated — not extracted)

# Lepetodrilus pustulosus — vent limpet Pacific; check for Atlantic Lepetodrilus
lepeto_ids <- rownames(tax_cur)[!is.na(tax_cur$Family) &
                                  tax_cur$Family == "Lepetodrilidae"]
write_fasta(lepeto_ids,
            "sequences_to_verify/Lepetodrilus_verify.fasta", tax_cur, ref_all)

# Strubellia / Naytia / Montfortula — Indo-Pacific genera; degraded to family
# Extract for record; re-verify if needed
strubellia_ids <- rownames(tax_cur)[!is.na(tax_cur$Phylum) &
                                      tax_cur$Phylum == "Mollusca" &
                                      !is.na(tax_cur$Family) &
                                      tax_cur$Family %in% c("Acochlidiidae",
                                                            "Nassariidae",
                                                            "Fissurellidae") &
                                      is.na(tax_cur$Genus)]
write_fasta(strubellia_ids,
            "sequences_to_verify/Strubellia_Naytia_Montfortula_verify.fasta", tax_cur, ref_all)

# Neritina — degraded to Gastropoda (class); verify if any Atlantic Neritina present
neritina_ids <- rownames(tax_cur)[!is.na(tax_cur$Family) &
                                    tax_cur$Family == "Neritidae" &
                                    is.na(tax_cur$Genus)]
write_fasta(neritina_ids,
            "sequences_to_verify/Neritina_verify.fasta", tax_cur, ref_all)

# ── 2. Cnidaria — DECIDED: Catostylus + Paraminabea eliminated (Step 25) ──────
# No sequences to extract.

# ── 3. Rhodophyta — DECIDED: all kept, degraded to Order level ────────────────
# Extract for record / future re-verification
rhodo_verify_ids <- rownames(tax_cur)[
  !is.na(tax_cur$Phylum) & tax_cur$Phylum == "Rhodophyta" &
    !is.na(tax_cur$Order) &
    tax_cur$Order %in% c("Gigartinales", "Rhodymeniales", "Corallinales",
                         "Gelidiales") &
    is.na(tax_cur$Genus)
]
write_fasta(rhodo_verify_ids,
            "sequences_to_verify/Rhodophyta_degraded_record.fasta", tax_cur, ref_all)

# ── 4. Ochrophyta — DECIDED: degraded to Fucales ─────────────────────────────
# Extract for record
ochro_verify_ids <- rownames(tax_cur)[
  !is.na(tax_cur$Phylum) & tax_cur$Phylum == "Ochrophyta" &
    !is.na(tax_cur$Order) & tax_cur$Order == "Fucales" &
    is.na(tax_cur$Genus)
]
write_fasta(ochro_verify_ids,
            "sequences_to_verify/Ochrophyta_Fucales_degraded_record.fasta", tax_cur, ref_all)

# ── 5. Annelida — Lumbrineridae reclassified OTU (for record) ─────────────────
lumbrinerid_otu <- "86102227f65ad9050ca97a771a11095e"
write_fasta(lumbrinerid_otu,
            "sequences_to_verify/Lumbrineridae_reclassified_record.fasta", tax_cur, ref_all)

# ── 6. Arthropoda — Gammaridae: 2 kept OTUs (for record) ─────────────────────
gammarus_keep_ids <- c("a98da7992e83dd20bbfcac8178d16534",
                       "bd0bf6019297e578ecc7fc44d8460a0e")
write_fasta(gammarus_keep_ids,
            "sequences_to_verify/Gammarus_kept_Gammaridae_record.fasta", tax_cur, ref_all)

cat("\nFASTA files written to: sequences_to_verify/\n")
cat("Files ending in '_verify' need further action.\n")
cat("Files ending in '_record' document decisions already applied.\n")


# ══════════════════════════════════════════════════════════════════════════════
# BARPLOT — Phylum composition (% reads)
# Faceted by: Locality × Host species × Tissue (KK / E)
#
# INPUT:  phyloseq_OTU97_LULU_e20_EUK_KKE_curated_speciesGlom.rds
# OUTPUT: barplot_phylum_locality_host_tissue.pdf  (landscape, ~A3)
#         barplot_phylum_locality_host_tissue.png  (300 dpi)
# ══════════════════════════════════════════════════════════════════════════════

library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(stringr)

setwd(".")  # set to your working directory containing phyloseq_OTU97_LULU_e20.rds

# ── 1. LOAD ───────────────────────────────────────────────────────────────────
ps <- readRDS("phyloseq_OTU97_LULU_e20_EUK_KKE_curated_speciesGlom.rds")
cat("Loaded:", ntaxa(ps), "taxa |", nsamples(ps), "samples\n")
cat("Sample variables:", paste(sample_variables(ps), collapse = ", "), "\n")

# ── 2. INFER HOST SPECIES from sample names if no dedicated column ─────────────
sd <- as.data.frame(sample_data(ps))

# Look for an existing host/species column
host_candidates <- c("host", "host_species", "Host", "Species", "urchin",
                     "host_sp", "Erizo", "HostSpecies", "especie", "ESPECIE")
host_col <- host_candidates[host_candidates %in% colnames(sd)]

if (length(host_col) == 0) {
  # Derive from sample name prefix:
  # AC* or AV* or AT* = Arbacia lixula
  # PC* or PV* or PT* = Paracentrotus lividus
  cat("  No host column found — inferring from sample name prefix.\n")
  sd$host_species <- case_when(
    str_detect(rownames(sd), "^A")  ~ "Arbacia lixula",
    str_detect(rownames(sd), "^P")  ~ "Paracentrotus lividus",
    TRUE                            ~ "Unknown"
  )
  sample_data(ps) <- sample_data(sd)
  host_col <- "host_species"
} else {
  host_col <- host_col[1]
}
cat("  Using host column:", host_col, "\n")
cat("  Host levels:", paste(unique(sd[[host_col]]), collapse = ", "), "\n")

# ── 3. COLLAPSE TO PHYLUM ─────────────────────────────────────────────────────
ps_phy <- tax_glom(ps, taxrank = "Phylum", NArm = TRUE)

# ── 4. MELT & COMPUTE RELATIVE ABUNDANCE PER SAMPLE ──────────────────────────
df_melt <- psmelt(ps_phy)

# rename host column to 'host_sp' using base R (avoids tidyselect issues)
colnames(df_melt)[colnames(df_melt) == host_col] <- "host_sp"

df_melt <- df_melt %>%
  mutate(
    tissue   = as.character(tissue),
    locality = as.character(locality),
    host_sp  = as.character(host_sp),
    Phylum   = as.character(Phylum)
  ) %>%
  # relative abundance within each sample
  group_by(Sample) %>%
  mutate(RelAbund = Abundance / sum(Abundance)) %>%
  ungroup()

# ── 5. MEAN RELATIVE ABUNDANCE PER GROUP ─────────────────────────────────────
df_group <- df_melt %>%
  group_by(locality, host_sp, tissue, Phylum) %>%
  summarise(
    MeanRelAbund = mean(RelAbund, na.rm = TRUE),
    n_samples    = n_distinct(Sample),
    .groups = "drop"
  ) %>%
  # re-normalise within group so bars sum to 1
  group_by(locality, host_sp, tissue) %>%
  mutate(MeanRelAbund = MeanRelAbund / sum(MeanRelAbund)) %>%
  ungroup()

# ── 6. COLOUR PALETTE — top N phyla + "Other" ─────────────────────────────────
N_TOP <- 15   # show top 15 phyla by total mean abundance; rest = "Other"

top_phyla <- df_group %>%
  group_by(Phylum) %>%
  summarise(total = sum(MeanRelAbund), .groups = "drop") %>%
  slice_max(total, n = N_TOP) %>%
  pull(Phylum)

df_plot <- df_group %>%
  mutate(Phylum_plot = ifelse(Phylum %in% top_phyla, Phylum, "Other")) %>%
  group_by(locality, host_sp, tissue, Phylum_plot) %>%
  summarise(MeanRelAbund = sum(MeanRelAbund), .groups = "drop") %>%
  group_by(locality, host_sp, tissue) %>%
  mutate(MeanRelAbund = MeanRelAbund / sum(MeanRelAbund)) %>%
  ungroup()

# Order phyla: most abundant first, "Other" last
phylum_order <- df_plot %>%
  filter(Phylum_plot != "Other") %>%
  group_by(Phylum_plot) %>%
  summarise(tot = sum(MeanRelAbund), .groups = "drop") %>%
  arrange(desc(tot)) %>%
  pull(Phylum_plot)
phylum_order <- c(phylum_order, "Other")

df_plot <- df_plot %>%
  mutate(Phylum_plot = factor(Phylum_plot, levels = rev(phylum_order)))

# Colour palette (up to 16 colours)
palette_base <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#D4A6C8", "#86BCB6", "#A0CBE8", "#FFBE7D", "#8CD17D",
  "#C7C7C7"   # Other
)
pal <- setNames(palette_base[seq_along(phylum_order)], phylum_order)

# ── 7. FACET LABELS ───────────────────────────────────────────────────────────
# tissue labels
tissue_labels <- c(KK = "Faeces (KK)", E = "Stomach (E)")

# Nice locality order (control sites first, vent last — adjust if needed)
locality_order <- sort(unique(df_plot$locality))
# If you know the order: locality_order <- c("Blanes","La Bajita","El Faro",...)

# host labels
host_labels <- c(
  "Arbacia lixula"       = "A. lixula",
  "Paracentrotus lividus"= "P. lividus"
)

df_plot <- df_plot %>%
  mutate(
    locality  = factor(locality, levels = locality_order),
    tissue    = factor(tissue,   levels = c("KK", "E"),
                       labels = c("Faeces (KK)", "Stomach (E)")),
    host_sp   = recode(host_sp, !!!host_labels),
    host_sp   = factor(host_sp)
  )

# ── 8. PLOT ───────────────────────────────────────────────────────────────────
p <- ggplot(df_plot,
            aes(x = host_sp, y = MeanRelAbund, fill = Phylum_plot)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.15) +
  facet_grid(tissue ~ locality, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = pal,
                    guide  = guide_legend(reverse = TRUE,
                                          ncol = 1,
                                          keywidth  = unit(0.4, "cm"),
                                          keyheight = unit(0.4, "cm"))) +
  labs(
    x    = NULL,
    y    = "Mean relative abundance (%)",
    fill = "Phylum",
    title = "Phylum composition by locality, host species and tissue type"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    strip.text         = element_text(face = "bold", size = 9),
    axis.text.x        = element_text(angle = 30, hjust = 1,
                                      face = "italic", size = 9),
    axis.text.y        = element_text(size = 8),
    axis.title.y       = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.title       = element_text(face = "bold", size = 9),
    legend.text        = element_text(size = 8),
    legend.position    = "right",
    plot.title         = element_text(size = 11, face = "bold",
                                      margin = margin(b = 8)),
    plot.margin        = margin(8, 8, 8, 8)
  )

# ── 9. SAVE ───────────────────────────────────────────────────────────────────
# Width scales with number of localities × 2 species
n_loc <- length(unique(df_plot$locality))
w_cm  <- max(20, n_loc * 6)

ggsave("barplot_phylum_locality_host_tissue.pdf",
       plot   = p,
       width  = w_cm,
       height = 18,
       units  = "cm",
       device = "pdf")

ggsave("barplot_phylum_locality_host_tissue.png",
       plot   = p,
       width  = w_cm,
       height = 18,
       units  = "cm",
       dpi    = 300)

cat("\nSaved:\n")
cat("  barplot_phylum_locality_host_tissue.pdf\n")
cat("  barplot_phylum_locality_host_tissue.png\n")

# ── 10. OPTIONAL: print summary table ─────────────────────────────────────────
cat("\nMean relative abundance (%) per group × phylum:\n")
df_group %>%
  filter(Phylum %in% top_phyla) %>%
  mutate(pct = round(MeanRelAbund * 100, 1)) %>%
  select(locality, host_sp, tissue, Phylum, pct, n_samples) %>%
  arrange(locality, host_sp, tissue, desc(pct)) %>%
  print(n = 80)