### ══════════════════════════════════════════════════════════════════════════
### OBJECTIVE 1 — Final analysis script (Blanes, faecal vs. stomach content)
###
### Faecal (F) vs Stomach content (S) — Blanes — Eukaryota only
### Dataset: phyloseq_OTU97_LULU_e20_EUK_KKE_curated_speciesGlom.rds
### (already collapsed at species level and cleaned of non-marine contaminants)
### Colours: Faecal (F) = orange (#F28E2B) | Stomach (S) = green (#59A14F)

# ── LIBRARIES ─────────────────────────────────────────────────────────────────
library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)
library(scales)
library(Matrix)
library(lme4)
library(lmerTest)
library(patchwork)
library(tibble)
library(permute)

# ── COLOUR PALETTE (fixed throughout) ────────────────────────────────────────
COL_F      <- "#F28E2B"   # orange — Faecal (F)
COL_S      <- "#59A14F"   # green  — Stomach (S)
COL_SHARED <- "#AAAAAA"   # grey   — shared

# ── LOAD OBJECT ───────────────────────────────────────────────────────────────
setwd(".")  # set to your working directory containing the final curated phyloseq object
ps <- readRDS("phyloseq_OTU97_LULU_e20_EUK_KKE_curated_speciesGlom.rds")

# ── COLUMN NAMES ──────────────────────────────────────────────────────────────
COL_LOCALITY   <- "locality"
COL_TISSUE     <- "tissue"       # values: "KK" (faecal) and "E" (stomach)
COL_SPECIES    <- "host_species"
COL_INDIVIDUAL <- "individual"
LOCALITY_NAME  <- "Blanes"

# ════════════════════════════════════════════════════════════════════════════════
# STEP 0 — Recode tissue labels: KK -> F, E -> S
# ════════════════════════════════════════════════════════════════════════════════
sample_data(ps)$sample_type <- dplyr::recode(
  as.character(sample_data(ps)[[COL_TISSUE]]),
  "KK" = "F",
  "E"  = "S"
)
COL_TYPE <- "sample_type"

# ════════════════════════════════════════════════════════════════════════════════
# STEP 1 — Filter Blanes + F and S samples
# ════════════════════════════════════════════════════════════════════════════════
ps_blanes <- subset_samples(ps, locality == LOCALITY_NAME)
ps_blanes <- prune_taxa(taxa_sums(ps_blanes) > 0, ps_blanes)

ps_blanes_FS <- subset_samples(ps_blanes, sample_type %in% c("F", "S"))
ps_blanes_FS <- prune_taxa(taxa_sums(ps_blanes_FS) > 0, ps_blanes_FS)

# ── PRE-RAREFACTION SUMMARY (for manuscript) ──────────────────────────────────
cat("\n=== PRE-RAREFACTION SUMMARY (Blanes, Eukaryota with host) ===\n")
cat("Total samples:", nsamples(ps_blanes_FS), "\n")
cat("  Faecal (F):  ", sum(sample_data(ps_blanes_FS)[[COL_TYPE]] == "F"), "\n")
cat("  Stomach (S): ", sum(sample_data(ps_blanes_FS)[[COL_TYPE]] == "S"), "\n")
cat("Total species (taxa):", ntaxa(ps_blanes_FS), "\n")
cat("Total eukaryotic reads:", sum(sample_sums(ps_blanes_FS)), "\n")

sd_pre  <- as.data.frame(sample_data(ps_blanes_FS))
otu_pre <- as(otu_table(ps_blanes_FS), "matrix")
if (!taxa_are_rows(ps_blanes_FS)) otu_pre <- t(otu_pre)
tax_pre <- as.data.frame(tax_table(ps_blanes_FS))

# Reads per host species and sample type
reads_summary <- data.frame(
  sample_type  = sd_pre[[COL_TYPE]],
  host_species = sd_pre[[COL_SPECIES]],
  reads        = colSums(otu_pre)
) %>%
  group_by(host_species, sample_type) %>%
  summarise(total_reads = sum(reads), n_samples = n(), .groups = "drop")
cat("\nReads per host species and sample type:\n")
print(reads_summary)

# Host DNA percentage per sample type
is_host <- grepl("Arbacia|Paracentrotus", tax_pre$Species, ignore.case = TRUE)
host_summary <- data.frame(
  sample_type = sd_pre[[COL_TYPE]],
  total       = colSums(otu_pre),
  host        = colSums(otu_pre[is_host, , drop = FALSE])
) %>%
  group_by(sample_type) %>%
  summarise(
    total_reads = sum(total),
    host_reads  = sum(host),
    pct_host    = round(sum(host) / sum(total) * 100, 1),
    .groups = "drop"
  )
cat("\nHost DNA percentage per sample type:\n")
print(host_summary)

cat("\nReads per sample (sorted):\n")
print(sort(sample_sums(ps_blanes_FS)))
cat("Minimum:", min(sample_sums(ps_blanes_FS)), "\n")
cat("Median: ", median(sample_sums(ps_blanes_FS)), "\n")

# ════════════════════════════════════════════════════════════════════════════════
# STEP 2 — Rarefaction curve
# ════════════════════════════════════════════════════════════════════════════════
otu_t <- as(otu_table(ps_blanes_FS), "matrix")
if (taxa_are_rows(ps_blanes_FS)) otu_t <- t(otu_t)

colores_rare <- ifelse(sample_data(ps_blanes_FS)[[COL_TYPE]] == "F", COL_F, COL_S)

pdf("OBJ1_EUK_rarefaction_curve_Blanes.pdf", width = 10, height = 6)
rarecurve(
  otu_t,
  step   = 200,
  sample = 90000,
  col    = colores_rare,
  label  = TRUE,
  cex    = 0.5,
  main   = "Rarefaction curve — Blanes (Eukaryota, species-level)"
)
abline(v = 90000, lty = 2, col = "grey40", lwd = 1.5)
legend("bottomright",
       legend = c("Faecal (F)", "Stomach (S)"),
       col    = c(COL_F, COL_S),
       lty = 1, lwd = 2, bty = "n")
dev.off()
cat("Saved: OBJ1_EUK_rarefaction_curve_Blanes.pdf\n")

# ════════════════════════════════════════════════════════════════════════════════
# STEP 3 — Rarefaction + Presence/Absence
# ════════════════════════════════════════════════════════════════════════════════
RARE_DEPTH <- 90000

set.seed(42)
ps_rare <- rarefy_even_depth(ps_blanes_FS,
                             sample.size = RARE_DEPTH,
                             rngseed     = 42,
                             replace     = FALSE)
ps_rare <- prune_taxa(taxa_sums(ps_rare) > 0, ps_rare)

cat("\n=== POST-RAREFACTION SUMMARY ===\n")
cat("Samples retained:", nsamples(ps_rare), "\n")
cat("Species retained:", ntaxa(ps_rare), "\n")
cat("Total reads:     ", sum(sample_sums(ps_rare)), "\n")
print(table(sample_data(ps_rare)[[COL_TYPE]]))

# Presence/Absence for beta diversity
ps_PA <- transform_sample_counts(ps_rare, function(x) ifelse(x > 0, 1, 0))
cat("Presence/Absence object ready\n")

# ════════════════════════════════════════════════════════════════════════════════
# STEP 4 — ALPHA DIVERSITY (Richness only)
# ════════════════════════════════════════════════════════════════════════════════
rich <- estimate_richness(ps_rare, measures = "Observed")
rich$SampleID     <- rownames(rich)
rich$sample_type  <- sample_data(ps_rare)[[COL_TYPE]]
rich$host_species <- sample_data(ps_rare)[[COL_SPECIES]]
rich$individual   <- sample_data(ps_rare)[[COL_INDIVIDUAL]]

cat("\n=== Richness by host species and sample type (mean +- SD) ===\n")
rich %>%
  group_by(host_species, sample_type) %>%
  summarise(mean = round(mean(Observed), 1),
            sd   = round(sd(Observed), 1),
            n    = n(), .groups = "drop") %>%
  print()

# ── Alpha boxplot ─────────────────────────────────────────────────────────────
p_alpha <- ggplot(rich, aes(x = sample_type, y = Observed, fill = sample_type)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8, width = 0.5) +
  geom_jitter(width = 0.12, alpha = 0.6, size = 2) +
  facet_wrap(~ host_species, labeller = label_wrap_gen(width = 25)) +
  scale_fill_manual(values = c("F" = COL_F, "S" = COL_S),
                    labels = c("F" = "Faecal (F)", "S" = "Stomach (S)")) +
  scale_x_discrete(labels = c("F" = "Faecal (F)", "S" = "Stomach (S)")) +
  labs(x    = "Sample type",
       y    = "Observed richness (no. species)",
       fill = "Sample type") +
  theme_bw(base_size = 13) +
  theme(strip.text       = element_text(face = "italic"),
        panel.grid.minor = element_blank(),
        legend.position  = "none")

ggsave("OBJ1_EUK_alpha_richness_Blanes_F_vs_S.pdf", p_alpha, width = 8, height = 5)
ggsave("OBJ1_EUK_alpha_richness_Blanes_F_vs_S.png", p_alpha, width = 8, height = 5, dpi = 300)
cat("Saved: alpha diversity plot\n")

# ── LMM full model ────────────────────────────────────────────────────────────
cat("\n=== LMM full model (sample_type * host_species | individual) ===\n")
lmm_full <- lmer(Observed ~ sample_type * host_species + (1 | individual), data = rich)
print(anova(lmm_full))

# ── LMM per species ───────────────────────────────────────────────────────────
for (sp in unique(rich$host_species)) {
  cat("\n--- LMM:", sp, "---\n")
  sub <- rich %>% filter(host_species == sp)
  m   <- lmer(Observed ~ sample_type + (1 | individual), data = sub)
  print(summary(m))
}

# ════════════════════════════════════════════════════════════════════════════════
# STEP 5 — BETA DIVERSITY: Jaccard PERMANOVA (constrained within individual)
#
# Differences in diet composition (beta diversity) between faecal and stomach
# content samples were assessed separately for each host species using Jaccard
# dissimilarity matrices and tested with PERMANOVA using the adonis2 function
# in the vegan package (Oksanen et al., 2022). Permutations were constrained
# within individual sea urchins to account for paired sampling.
# ════════════════════════════════════════════════════════════════════════════════
jac_dist <- phyloseq::distance(ps_PA, method = "jaccard", binary = TRUE)
ord_jac  <- ordinate(ps_PA, method = "PCoA", distance = jac_dist)

sd_df    <- as.data.frame(sample_data(ps_PA))
sd_clean <- data.frame(
  sample_type  = as.factor(sd_df[[COL_TYPE]]),
  host_species = as.factor(sd_df[[COL_SPECIES]]),
  individual   = as.factor(sd_df[[COL_INDIVIDUAL]]),
  row.names    = rownames(sd_df)
)
jac_dist2 <- as.dist(as.matrix(jac_dist)[rownames(sd_clean), rownames(sd_clean)])

# ── PERMANOVA per species ─────────────────────────────────────────────────────
run_permanova_species <- function(species_name, dist_full, sd_full) {
  cat("\n=== PERMANOVA Jaccard —", species_name, "===\n")
  idx     <- sd_full$host_species == species_name
  samps   <- rownames(sd_full)[idx]
  sd_sp   <- droplevels(sd_full[idx, ])
  dist_sp <- as.dist(as.matrix(dist_full)[samps, samps])
  perm_ctrl <- permute::how(blocks = sd_sp$individual, nperm = 999)
  res <- adonis2(dist_sp ~ sample_type, data = sd_sp, permutations = perm_ctrl)
  print(res)
  return(res)
}

perm_arb <- run_permanova_species("Arbacia lixula",        jac_dist2, sd_clean)
perm_par <- run_permanova_species("Paracentrotus lividus", jac_dist2, sd_clean)

# ── Betadisper per species ────────────────────────────────────────────────────
run_betadisper_species <- function(species_name, dist_full, sd_full) {
  cat("\n=== Betadisper —", species_name, "===\n")
  idx     <- sd_full$host_species == species_name
  samps   <- rownames(sd_full)[idx]
  sd_sp   <- droplevels(sd_full[idx, ])
  dist_sp <- as.dist(as.matrix(dist_full)[samps, samps])
  bd      <- betadisper(dist_sp, sd_sp$sample_type)
  res     <- permutest(bd)
  print(res)
  cat("Group distances:\n")
  print(bd$group.distances)
  return(res)
}

bd_arb <- run_betadisper_species("Arbacia lixula",        jac_dist2, sd_clean)
bd_par <- run_betadisper_species("Paracentrotus lividus", jac_dist2, sd_clean)

# ── PCoA plot per species ─────────────────────────────────────────────────────
plot_pcoa_species <- function(species_name, dist_full, sd_full, perm_res, bd_res) {
  
  idx     <- sd_full$host_species == species_name
  samps   <- rownames(sd_full)[idx]
  sd_sp   <- droplevels(sd_full[idx, ])
  dist_sp <- as.dist(as.matrix(dist_full)[samps, samps])
  
  ord_sp  <- cmdscale(dist_sp, eig = TRUE, k = 2)
  eig     <- ord_sp$eig
  pct     <- round(eig / sum(eig[eig > 0]) * 100, 1)
  
  pcoa_df <- as.data.frame(ord_sp$points)
  colnames(pcoa_df) <- c("Axis1", "Axis2")
  pcoa_df$sample_type <- sd_sp[rownames(pcoa_df), "sample_type"]
  pcoa_df$individual  <- sd_sp[rownames(pcoa_df), "individual"]
  
  sp_label   <- ifelse(grepl("Arbacia", species_name), "A. lixula", "P. lividus")
  perm_p     <- perm_res$`Pr(>F)`[1]
  perm_p_txt <- ifelse(perm_p < 0.05, paste0("p = ", round(perm_p, 3)),
                       paste0("p = ", round(perm_p, 3), " (n.s.)"))
  bd_p       <- bd_res$tab$`Pr(>F)`[1]
  bd_p_txt   <- ifelse(bd_p < 0.05, paste0("p = ", round(bd_p, 3)),
                       paste0("p = ", round(bd_p, 3), " (n.s.)"))
  
  perm_label <- paste0(
    "PERMANOVA: R\u00b2 = ", round(perm_res$R2[1], 3), ", ", perm_p_txt, "\n",
    "Betadisper: F = ",      round(bd_res$tab$F[1], 2), ", ", bd_p_txt
  )
  
  ggplot(pcoa_df, aes(x = Axis1, y = Axis2)) +
    stat_ellipse(aes(color = sample_type, group = sample_type),
                 linetype = 2, linewidth = 0.8, level = 0.95) +
    geom_point(aes(color = sample_type), size = 4, alpha = 0.9) +
    geom_line(aes(group = individual),
              color = "grey60", linewidth = 0.4, alpha = 0.6) +
    scale_color_manual(values = c("F" = COL_F, "S" = COL_S),
                       labels = c("F" = "Faecal (F)", "S" = "Stomach (S)")) +
    annotate("text", x = Inf, y = Inf, label = perm_label,
             hjust = 1.05, vjust = 1.5, size = 3, color = "grey30") +
    labs(x     = paste0("PCoA Axis 1 (", pct[1], "%)"),
         y     = paste0("PCoA Axis 2 (", pct[2], "%)"),
         title = bquote(italic(.(sp_label))),
         color = "Sample type") +
    theme_bw(base_size = 13) +
    theme(plot.title       = element_text(hjust = 0.5, face = "bold.italic"),
          panel.grid.minor = element_blank())
}

p_pcoa_arb <- plot_pcoa_species("Arbacia lixula",        jac_dist2, sd_clean, perm_arb, bd_arb)
p_pcoa_par <- plot_pcoa_species("Paracentrotus lividus", jac_dist2, sd_clean, perm_par, bd_par)

p_pcoa_combined <- p_pcoa_arb + p_pcoa_par +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 14, face = "bold"))
  )

ggsave("OBJ1_EUK_beta_jaccard_PCoA_by_species.pdf", p_pcoa_combined, width = 12, height = 5)
ggsave("OBJ1_EUK_beta_jaccard_PCoA_by_species.png", p_pcoa_combined, width = 12, height = 5, dpi = 300)
cat("Saved: PCoA per species\n")

# ════════════════════════════════════════════════════════════════════════════════
# STEP 6 — VENN DIAGRAM: Faecal (F) vs Stomach (S) per species
# ════════════════════════════════════════════════════════════════════════════════
otu_pa <- as(otu_table(ps_PA), "matrix")
if (!taxa_are_rows(ps_PA)) otu_pa <- t(otu_pa)
sd_pa  <- as.data.frame(sample_data(ps_PA))

get_taxa <- function(type_val, species_val) {
  samps <- rownames(sd_pa)[sd_pa[[COL_TYPE]] == type_val &
                             sd_pa[[COL_SPECIES]] == species_val]
  rownames(otu_pa)[rowSums(otu_pa[, samps, drop = FALSE]) > 0]
}

plot_venn_ggplot <- function(species_name, filename) {
  
  taxa_F <- get_taxa("F", species_name)
  taxa_S <- get_taxa("S", species_name)
  
  only_F <- length(setdiff(taxa_F, taxa_S))
  only_S <- length(setdiff(taxa_S, taxa_F))
  shared <- length(intersect(taxa_F, taxa_S))
  pct_F  <- round(shared / length(taxa_F) * 100, 1)
  pct_S  <- round(shared / length(taxa_S) * 100, 1)
  
  cat("\n===", species_name, "===\n")
  cat("Only Faecal (F): ", only_F,  "\n")
  cat("Only Stomach (S):", only_S,  "\n")
  cat("Shared:          ", shared,  "\n")
  cat("% of F species shared with S:  ", pct_F, "%\n")
  cat("% of S species recovered in F: ", pct_S, "%\n")
  
  cx_F <- -0.28; cx_S <- 0.28; cy <- 0.00; r_circ <- 0.52
  theta  <- seq(0, 2 * pi, length.out = 300)
  circ_F <- data.frame(x = cx_F + r_circ * cos(theta),
                       y = cy   + r_circ * sin(theta))
  circ_S <- data.frame(x = cx_S + r_circ * cos(theta),
                       y = cy   + r_circ * sin(theta))
  
  label_df <- data.frame(
    x     = c(-0.60, 0.00, 0.60),
    y     = c( 0.00, 0.00, 0.00),
    label = c(as.character(only_F), as.character(shared), as.character(only_S))
  )
  set_df <- data.frame(
    x = c(cx_F, cx_S), y = c(0.68, 0.68),
    label = c("Faecal (F)", "Stomach (S)")
  )
  sub_txt <- paste0(
    shared, " shared species  |  ",
    pct_S, "% of stomach species recovered in faeces  |  ",
    pct_F, "% of faecal species shared with stomach"
  )
  sp_label <- ifelse(grepl("Arbacia", species_name), "A. lixula", "P. lividus")
  
  p <- ggplot() +
    geom_polygon(data = circ_F, aes(x = x, y = y),
                 fill = COL_F, alpha = 0.35, color = COL_F, linewidth = 0.8) +
    geom_polygon(data = circ_S, aes(x = x, y = y),
                 fill = COL_S, alpha = 0.35, color = COL_S, linewidth = 0.8) +
    geom_text(data = label_df, aes(x = x, y = y, label = label),
              size = 8, fontface = "bold", color = "grey20") +
    geom_text(data = set_df, aes(x = x, y = y, label = label),
              size = 4, fontface = "bold", color = c(COL_F, COL_S)) +
    labs(title = bquote(italic(.(sp_label))), subtitle = sub_txt) +
    coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-0.85, 0.85)) +
    theme_void(base_size = 12) +
    theme(
      plot.title    = element_text(hjust = 0.5, size = 13, face = "bold.italic",
                                   margin = margin(b = 4)),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey40",
                                   margin = margin(b = 6)),
      plot.margin   = margin(10, 10, 10, 10)
    )
  
  ggsave(paste0(filename, ".pdf"), p, width = 5, height = 5)
  ggsave(paste0(filename, ".png"), p, width = 5, height = 5, dpi = 300)
  cat("Saved:", filename, "\n")
  return(p)
}

p_arb <- plot_venn_ggplot("Arbacia lixula",       "OBJ1_venn_Arbacia_F_vs_S")
p_par <- plot_venn_ggplot("Paracentrotus lividus", "OBJ1_venn_Paracentrotus_F_vs_S")

p_venn_combined <- p_arb + p_par +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag   = element_text(size = 14, face = "bold"))
  )

ggsave("OBJ1_venn_combined.pdf", p_venn_combined, width = 14, height = 5)
ggsave("OBJ1_venn_combined.png", p_venn_combined, width = 14, height = 6, dpi = 300)
cat("Saved: combined Venn panel\n")

# ════════════════════════════════════════════════════════════════════════════════
# STEP 7 — Species exclusive to Stomach (S): taxonomy table
# ════════════════════════════════════════════════════════════════════════════════
get_exclusive_S <- function(species_name) {
  
  taxa_F      <- get_taxa("F", species_name)
  taxa_S      <- get_taxa("S", species_name)
  exclusive_S <- setdiff(taxa_S, taxa_F)
  
  tax      <- as.data.frame(tax_table(ps_PA))
  tax_excl <- tax[exclusive_S, , drop = FALSE] %>% rownames_to_column("OTU_ID")
  
  otu_mat <- as(otu_table(ps_PA), "matrix")
  if (!taxa_are_rows(ps_PA)) otu_mat <- t(otu_mat)
  
  sd_pa2 <- as.data.frame(sample_data(ps_PA))
  samp_S <- rownames(sd_pa2)[sd_pa2[[COL_TYPE]] == "S" &
                               sd_pa2[[COL_SPECIES]] == species_name]
  
  reads_S <- otu_mat[exclusive_S, samp_S, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("OTU_ID") %>%
    mutate(
      total_reads_S     = rowSums(across(where(is.numeric))),
      n_samples_present = rowSums(across(where(is.numeric)) > 0)
    )
  
  result <- tax_excl %>%
    left_join(reads_S %>% select(OTU_ID, total_reads_S, n_samples_present),
              by = "OTU_ID") %>%
    arrange(desc(total_reads_S))
  
  cat("\n===", species_name, "— species exclusive to Stomach (S) ===\n")
  print(result %>% select(OTU_ID, Phylum, Class, Order, Family, Genus, Species,
                          total_reads_S, n_samples_present))
  return(result)
}

excl_S_arbacia <- get_exclusive_S("Arbacia lixula")
excl_S_paracen <- get_exclusive_S("Paracentrotus lividus")

write.csv(excl_S_arbacia, "OBJ1_EUK_exclusive_S_Arbacia.csv",       row.names = FALSE)
write.csv(excl_S_paracen, "OBJ1_EUK_exclusive_S_Paracentrotus.csv", row.names = FALSE)
cat("Saved: exclusive stomach species tables\n")

# ════════════════════════════════════════════════════════════════════════════════
# STEP 8 — Species exclusive to Faecal (F): taxonomy table
# ════════════════════════════════════════════════════════════════════════════════
get_exclusive_F <- function(species_name) {
  
  taxa_F      <- get_taxa("F", species_name)
  taxa_S      <- get_taxa("S", species_name)
  exclusive_F <- setdiff(taxa_F, taxa_S)
  
  tax      <- as.data.frame(tax_table(ps_PA))
  tax_excl <- tax[exclusive_F, , drop = FALSE] %>% rownames_to_column("OTU_ID")
  
  otu_mat <- as(otu_table(ps_rare), "matrix")   # use counts (not PA) for reads
  if (!taxa_are_rows(ps_rare)) otu_mat <- t(otu_mat)
  
  sd_pa2 <- as.data.frame(sample_data(ps_PA))
  samp_F <- rownames(sd_pa2)[sd_pa2[[COL_TYPE]] == "F" &
                               sd_pa2[[COL_SPECIES]] == species_name]
  
  reads_F <- otu_mat[exclusive_F, samp_F, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("OTU_ID") %>%
    mutate(
      total_reads_F     = rowSums(across(where(is.numeric))),
      n_samples_present = rowSums(across(where(is.numeric)) > 0)
    )
  
  result <- tax_excl %>%
    left_join(reads_F %>% select(OTU_ID, total_reads_F, n_samples_present),
              by = "OTU_ID") %>%
    arrange(desc(total_reads_F))
  
  cat("\n===", species_name, "— species exclusive to Faecal (F) ===\n")
  print(result %>% select(OTU_ID, Phylum, Class, Order, Family, Genus, Species,
                          total_reads_F, n_samples_present))
  return(result)
}

excl_F_arbacia <- get_exclusive_F("Arbacia lixula")
excl_F_paracen <- get_exclusive_F("Paracentrotus lividus")

write.csv(excl_F_arbacia, "OBJ1_EUK_exclusive_F_Arbacia.csv",       row.names = FALSE)
write.csv(excl_F_paracen, "OBJ1_EUK_exclusive_F_Paracentrotus.csv", row.names = FALSE)
cat("Saved: exclusive faecal species tables\n")

rich %>%
  group_by(host_species, sample_type) %>%
  summarise(mean = round(mean(Observed),1), sd = round(sd(Observed),1), n=n())

cat("\n✅ OBJECTIVE 1 COMPLETE — all files saved\n")
