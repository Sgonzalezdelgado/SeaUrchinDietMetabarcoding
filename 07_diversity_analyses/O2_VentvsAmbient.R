### ══════════════════════════════════════════════════════════════════════════
### OBJECTIVE 2 — Final analysis script (La Palma, vent vs ambient localities)
###   0. Setup
###   1. Load & filter (host + Echinodermata removed, Los Porretos excluded)
###   2. Rarefaction curve (min_reads = 5000)
###   3. Alpha diversity (LMM) -> Table S5, Fig. 3A
###   4. Beta diversity, PCoA panel -> Fig. 3C-D
###   4b. PERMANOVA by Condition/Locality -> Table S6A-C
###   5. PERMANOVA by Phylum + Family (FDR) -> Table S7
###   7. IndVal.g by locality (mixed-rank: species > genus > family > phylum,
###      whichever is most resolved per taxon) + dotplot coloured by
###      locality -> Fig. 4, Table S8
###
## ---- 0. SETUP -----------------------------------------------

setwd(".")  # set to your working directory containing the final curated phyloseq object

pkgs_cran <- c("vegan", "ggplot2", "dplyr", "tidyr", "stringr",
               "tibble", "scales", "patchwork", "lme4", "lmerTest",
               "emmeans", "car", "indicspecies",
               "RColorBrewer")
pkgs_bio  <- c("phyloseq", "Biostrings")

for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in pkgs_bio) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)
}
invisible(lapply(c(pkgs_cran, pkgs_bio), library, character.only = TRUE))

OUTDIR <- "Objetivo2_v3"
dir.create(OUTDIR, showWarnings = FALSE)
cat("Outputs en:", file.path(getwd(), OUTDIR), "\n")


## ---- 1. CARGA Y FILTROS -------------------------------------
cat("\n=== 1. Carga y filtros ===\n")

# Dataset nuevo: ya colapsado por especie, solo Eukariota, KK + E curados
ps_full <- readRDS("phyloseq_OTU97_LULU_e20_EUK_KKE_curated_speciesGlom.rds")
cat("Phyloseq cargado:\n"); print(ps_full)

# -- Detectar columnas de metadata -----------------------------------
sd_all <- as.data.frame(sample_data(ps_full))
cat("\nColumnas de metadata disponibles:\n"); print(colnames(sd_all))

kk_col   <- c("tissue","sample_type","matrix","Material","type")[
  c("tissue","sample_type","matrix","Material","type") %in% colnames(sd_all)][1]
cond_col <- c("condition","CONDICION","Condition","pH_level","treatment")[
  c("condition","CONDICION","Condition","pH_level","treatment") %in% colnames(sd_all)][1]
loc_col  <- c("locality","Localidad","location","Site","Locality")[
  c("locality","Localidad","location","Site","Locality") %in% colnames(sd_all)][1]
sp_col   <- c("host_species","ESPECIE","species","Erizo","Species")[
  c("host_species","ESPECIE","species","Erizo","Species") %in% colnames(sd_all)][1]

cat("\nColumnas detectadas — KK:", kk_col, "| Cond:", cond_col,
    "| Loc:", loc_col, "| Sp:", sp_col, "\n")
stopifnot(!is.na(c(kk_col, cond_col, loc_col, sp_col)))

# -- 1a) Solo KK -------------------------------------------------
ps <- subset_samples(ps_full, get(kk_col) == "KK")
cat("Tras filtrar KK:", nsamples(ps), "muestras\n")

# -- 1b) Solo La Palma (excluir Blanes) --------------------------
ps <- subset_samples(ps, !grepl("Blanes|blanes", get(loc_col), ignore.case = TRUE))
cat("Tras excluir Blanes:", nsamples(ps), "muestras\n")

# -- 1c) Armonizar nombres Transition -> Vent / Los Porretos -----
sd_fix <- as.data.frame(sample_data(ps))
sd_fix[[loc_col]]  <- gsub("Transition Faro|Transicion Faro|TransicionFaro",
                           "Los Porretos", sd_fix[[loc_col]])
sd_fix[[cond_col]] <- gsub("Transition|Transicion",
                           "Vent", sd_fix[[cond_col]])
sample_data(ps) <- sample_data(sd_fix)

# -- 1d) Eliminar Los Porretos (n insuficiente) ------------------
ps <- subset_samples(ps, !grepl("Los Porretos|LaPorretos", get(loc_col)))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
cat("Tras eliminar Los Porretos:", nsamples(ps), "muestras\n")

# -- 1e) Eliminar host (erizo) y todo Echinodermata -----------------------------------
tax_df <- as.data.frame(tax_table(ps))

host_idx <- grepl("Arbacia|Paracentrotus", tax_df$Species, ignore.case = TRUE) |
  (tax_df$Phylum == "Echinodermata" & !is.na(tax_df$Phylum))

if (sum(host_idx) > 0) {
  ps <- prune_taxa(!host_idx, ps)
  cat("OTUs eliminados (host + Echinodermata):", sum(host_idx), "\n")
} else {
  cat("Sin OTUs host/Echinodermata en dataset.\n")
}
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
cat("\n=== Objeto tras filtros iniciales ===\n"); print(ps)

# -- Distribución -------------------------------------------------
sd_ps <- as.data.frame(as.matrix(sample_data(ps)))
cat("\nDistribución Condición × Localidad × Especie:\n")
print(table(sd_ps[[loc_col]], sd_ps[[cond_col]], sd_ps[[sp_col]]))

saveRDS(ps, file.path(OUTDIR, "ps_KK_LaPalma_filtered.rds"))

# -- Separar por especie -----------------------------------------
ps_arb <- subset_samples(ps, get(sp_col) %in%
                           grep("Arbacia", unique(sd_ps[[sp_col]]),
                                value = TRUE, ignore.case = TRUE))
ps_arb <- prune_taxa(taxa_sums(ps_arb) > 0, ps_arb)

ps_par <- subset_samples(ps, get(sp_col) %in%
                           grep("Paracentrotus", unique(sd_ps[[sp_col]]),
                                value = TRUE, ignore.case = TRUE))
ps_par <- prune_taxa(taxa_sums(ps_par) > 0, ps_par)

cat("\nArbacia lixula:", nsamples(ps_arb), "muestras\n")
cat("Paracentrotus lividus:", nsamples(ps_par), "muestras\n")

ps_list  <- list("Arbacia_lixula"       = ps_arb,
                 "Paracentrotus_lividus" = ps_par)
sp_names <- list("Arbacia_lixula"       = "Arbacia lixula",
                 "Paracentrotus_lividus" = "Paracentrotus lividus")

## Paletas --------------------------------------------------------
pal_loc <- c(
  "Cercana La Bajita" = "#1B9E77",
  "Charco Verde"      = "#80B1D3",
  "El Faro"           = "#FDB462",
  "La Bajita"         = "#66C2A5",
  "La Bombilla"       = "#FB8072"
)

loc_order_all <- c("La Bajita", "Cercana La Bajita","Charco Verde",
                   "El Faro", "La Bombilla")

pal_phyla <- c(
  "Rhodophyta"      = "#C0392B",   # rojo
  "Ochrophyta"      = "#E67E22",   # naranja
  "Annelida"        = "#9B59B6",   # morado
  "Arthropoda"      = "#2980B9",   # azul oscuro
  "Cnidaria"        = "#5DADE2",   # azul claro
  "Nemertea"        = "#1E8449",   # verde oscuro
  "Porifera"        = "#F39C12",   # amarillo-naranja
  "Mollusca"        = "#1ABC9C",   # turquesa
  "Platyhelminthes" = "#2ECC71",   # verde claro
  "Chlorophyta"     = "#A9DFBF",   # verde menta
  "Echinodermata"   = "#F1948A",   # rosa salmón
  "Oomycota"        = "#D7BDE2",   # lila claro
  "Chordata"        = "#FAD7A0",   # melocotón
  "Bacillariophyta" = "#AED6F1",   # azul muy claro
  "Haptophyta"      = "#A3E4D7",   # verde agua claro
  "Bryozoa"         = "#CCD1D1",   # gris claro
  "Amoebozoa"       = "#797D7F"    # gris oscuro
)

## ---- 2. CURVA DE RAREFACCIÓN ----------------------------------
cat("\n=== 2. Curvas de rarefacción ===\n")

min_reads <- 5000

rarecurve_df <- function(otu_mat, step = 200) {
  bind_rows(lapply(colnames(otu_mat), function(s) {
    x  <- otu_mat[, s]; x <- x[x > 0]; n <- sum(x)
    ns <- unique(c(seq(0, n, by = step), n))
    data.frame(
      Sample = s, Reads = ns,
      OTUs = sapply(ns, function(ni) {
        if (ni == 0) return(0)
        tryCatch(vegan::rarefy(x, ni), error = function(e) NA_real_)
      })
    )
  }))
}

for (sp_label in names(ps_list)) {
  
  ps_sp   <- ps_list[[sp_label]]
  sp_name <- sp_names[[sp_label]]
  
  otu_sp <- as(otu_table(ps_sp), "matrix")
  if (!taxa_are_rows(ps_sp)) otu_sp <- t(otu_sp)
  
  cat("\n---", sp_name, "---\n")
  reads_sp <- sort(colSums(otu_sp))
  cat("Min:", min(reads_sp), "| Mediana:", median(reads_sp),
      "| < 5000:", sum(reads_sp < min_reads), "muestras excluidas\n")
  
  meta_sp <- as.data.frame(as.matrix(sample_data(ps_sp))) %>%
    rownames_to_column("Sample") %>%
    select(Sample, all_of(c(cond_col, loc_col)))
  
  rc <- rarecurve_df(otu_sp, step = 200) %>%
    left_join(meta_sp, by = "Sample")
  
  p_rc <- ggplot(rc, aes(x = Reads, y = OTUs,
                         group = Sample,
                         colour = .data[[cond_col]])) +
    geom_line(alpha = 0.45, linewidth = 0.4) +
    geom_vline(xintercept = min_reads, linetype = "dashed",
               linewidth = 0.8, colour = "black") +
    annotate("text", x = min_reads * 1.1,
             y = max(rc$OTUs, na.rm = TRUE) * 0.92,
             label = paste0(min_reads, " reads (rarefaction threshold)"),
             hjust = 0, size = 3.5) +
    scale_colour_manual(
      values = c("Control" = "#2166AC", "Vent" = "#D73027")) +
    scale_x_continuous(labels = scales::comma, limits = c(0, 20000)) +
    labs(title    = paste("Rarefaction curve –", sp_name),
         subtitle = "KK samples | La Palma",
         x = "Number of reads", y = "Observed OTUs",
         colour = "Condition") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "italic"))
  
  ggsave(file.path(OUTDIR, paste0("Fig_rarefaction_", sp_label, ".pdf")),
         plot = p_rc, width = 9, height = 6)
  ggsave(file.path(OUTDIR, paste0("Fig_rarefaction_", sp_label, ".png")),
         plot = p_rc, width = 9, height = 6, dpi = 300)
  cat("Guardado: Fig_rarefaction_", sp_label, "\n")
}

##---- 3. ALPHA DIVERSITY (LMM) ---------------------------------
## Riqueza ~ pH + (1|Locality)
## ---------------------------------------------------------------
cat("\n=== 3. Alpha diversity — LMM ===\n")
 
sink(file.path(OUTDIR, "Stats_alpha_LMM.txt"))
cat("=== ALPHA DIVERSITY: Observed OTU richness ===\n")
cat("Rarefaction to", min_reads, "reads\n")
cat("Model: Richness ~ pH_condition + (1|Locality)\n")
sink()
 
## ---- 3a. ALPHA DIVERSITY — ESTADÍSTICA (LMM) ----
 
run_alpha_stats <- function(sp_label) {
  
  ps_sp   <- ps_list[[sp_label]]
  sp_name <- sp_names[[sp_label]]
  cat("\n---", sp_name, "---\n")
  
  # Rarefacción; excluye muestras < min_reads
  cat("Muestras con <", min_reads, "reads (excluidas):",
      sum(sample_sums(ps_sp) < min_reads), "\n")
  set.seed(42)
  ps_r <- rarefy_even_depth(ps_sp, sample.size = min_reads,
                            rngseed = 42, replace = FALSE,
                            trimOTUs = TRUE, verbose = FALSE)
  cat("Muestras tras rarefacción:", nsamples(ps_r), "\n")
  
  # Construir alpha_df en dos pasos separados
  rich  <- estimate_richness(ps_r, measures = "Observed")
  sd_df <- as.data.frame(as.matrix(sample_data(ps_r))) %>%
    rownames_to_column("Sample")
  
  alpha_df <- rich %>%
    rownames_to_column("Sample") %>%
    left_join(sd_df, by = "Sample") %>%
    mutate(
      Observed    = as.numeric(Observed),
      Condition_f = as.factor(.data[[cond_col]]),
      Locality_f  = as.factor(.data[[loc_col]])
    )
  
  # Renombrar Control -> Ambient
  levels(alpha_df$Condition_f) <- gsub("Control", "Ambient",
                                       levels(alpha_df$Condition_f))
  
  cat("Niveles Condition_f:", levels(alpha_df$Condition_f), "\n")
  cat("Niveles Locality_f:",  levels(alpha_df$Locality_f),  "\n")
  
  # ---- LMM: Richness ~ Condition + (1|Locality) ----
  lmm <- lmer(Observed ~ Condition_f + (1|Locality_f),
              data = alpha_df, REML = TRUE)
  
  if (isSingular(lmm)) cat("  AVISO: modelo singular (varianza Locality ~ 0)\n")
  
  cat("\n--- LMM: Richness ~ Condition + (1|Locality) ---\n")
  print(summary(lmm))
  print(anova(lmm, ddf = "Kenward-Roger"))
  print(summary(pw_loc))
  
  # ---- Guardar stats ----
  sink(file.path(OUTDIR, "Stats_alpha_LMM.txt"), append = TRUE)
  cat("\n\n===", sp_name, "===\n")
  cat("--- LMM summary ---\n\n")
  print(summary(lmm))
  cat("\n--- ANOVA (Kenward-Roger) ---\n\n")
  print(anova(lmm, ddf = "Kenward-Roger"))
  sink()
  
  list(
    sp_label = sp_label, sp_name = sp_name,
    alpha_df = alpha_df, lmm = lmm
  )
}

 
## ---- 3b. ALPHA DIVERSITY — GRÁFICO (boxplot por localidad) ----
 
plot_alpha_boxplot <- function(res) {
  
  alpha_df <- res$alpha_df
  sp_label <- res$sp_label
  sp_name  <- res$sp_name
  
  loc_order_sp <- loc_order_all[loc_order_all %in% levels(alpha_df$Locality_f)]
  alpha_df$Locality_f <- factor(alpha_df$Locality_f, levels = loc_order_sp)
  n_ctrl_sp <- sum(loc_order_sp %in% c("La Bajita", "Charco Verde", "El Hoyo"))
  
  p_alpha <- ggplot(alpha_df,
                    aes(x = Locality_f, y = Observed, fill = Condition_f)) +
    annotate("rect", xmin = 0.5, xmax = n_ctrl_sp + 0.5,
             ymin = -Inf, ymax = Inf, fill = "#DDEEFF", alpha = 0.3) +
    annotate("rect", xmin = n_ctrl_sp + 0.5,
             xmax = length(loc_order_sp) + 0.5,
             ymin = -Inf, ymax = Inf, fill = "#FFDDDD", alpha = 0.3) +
    geom_boxplot(outlier.shape = 21, width = 0.6, linewidth = 0.4) +
    geom_jitter(aes(colour = Locality_f), width = 0.1,
                size = 2.2, alpha = 0.7) +
    geom_vline(xintercept = n_ctrl_sp + 0.5,
               linetype = "dashed", colour = "grey40", linewidth = 0.6) +
    annotate("text", x = n_ctrl_sp / 2 + 0.5, y = Inf,
             label = "AMBIENT", vjust = 1.5, size = 3.5,
             colour = "#2166AC", fontface = "bold") +
    annotate("text",
             x = n_ctrl_sp + (length(loc_order_sp) - n_ctrl_sp) / 2 + 0.5,
             y = Inf, label = "VENT", vjust = 1.5, size = 3.5,
             colour = "#D73027", fontface = "bold") +
    scale_fill_manual(values = c("Ambient" = "#7BAFD4",
                                 "Vent"    = "#F4A582")) +
    scale_colour_manual(values = pal_loc) +
    scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
    labs(title    = paste("Alpha diversity (Observed richness) \u2013", sp_name),
         subtitle = paste("Rarefied to", min_reads,
                          "reads | KK samples | La Palma"),
         x = "Locality", y = "Observed OTUs",
         fill = "Condition", colour = "Locality") +
    theme_bw(base_size = 13) +
    theme(legend.position  = "right",
          plot.title       = element_text(face = "italic"),
          axis.text.x      = element_text(size = 9))
  
  print(p_alpha)
  ggsave(file.path(OUTDIR, paste0("Fig_alpha_richness_", sp_label, ".pdf")),
         plot = p_alpha, width = 10, height = 6)
  ggsave(file.path(OUTDIR, paste0("Fig_alpha_richness_", sp_label, ".png")),
         plot = p_alpha, width = 10, height = 6, dpi = 300)
  cat("Guardado: Fig_alpha_richness_", sp_label, "\n")
  
  invisible(p_alpha)
}
 
res_arb  <- run_alpha_stats("Arbacia_lixula")
plot_alpha_boxplot(res_arb)
res_par  <- run_alpha_stats("Paracentrotus_lividus")
plot_alpha_boxplot(res_par)

## ---- 4. BETA DIVERSITY (Jaccard P/A) --------------------------
## PERMANOVA: adonis2 ~ Locality (locality as fixed factor)
## Pairwise entre localidades (FDR)
## Betadisper para homogeneidad de dispersión
## ---------------------------------------------------------------
cat("\n=== 4. Beta diversity — Jaccard PERMANOVA ===\n")
 
sink(file.path(OUTDIR, "Stats_beta_Jaccard.txt"))
cat("=== BETA DIVERSITY: Jaccard (presence/absence) ===\n")
cat("Model: adonis2(Jaccard ~ Locality)\n")
cat("Pairwise: between localities, FDR correction\n")
cat("Dispersion: betadisper by Condition and by Locality\n\n")
sink()
 
 
## ---- 4. BETA DIVERSITY — ESTADÍSTICA (PERMANOVA + pairwise + betadisper) ----
 

run_beta_stats <- function(sp_label) {
  
  ps_sp   <- ps_list[[sp_label]]
  sp_name <- sp_names[[sp_label]]
  cat("\n---", sp_name, "---\n")
  
  # Rarefy
  set.seed(42)
  ps_r <- rarefy_even_depth(ps_sp, sample.size = min_reads,
                            rngseed = 42, replace = FALSE,
                            trimOTUs = TRUE, verbose = FALSE)
  
  otu_r <- as(otu_table(ps_r), "matrix")
  if (!taxa_are_rows(ps_r)) otu_r <- t(otu_r)
  otu_t <- t(otu_r)
  
  jac <- vegdist(otu_t, method = "jaccard", binary = TRUE)
  
  sd_sp <- as.data.frame(as.matrix(sample_data(ps_r)))
  sd_sp[[cond_col]] <- as.factor(sd_sp[[cond_col]])
  sd_sp[[loc_col]]  <- as.factor(sd_sp[[loc_col]])
  
  # Renombrar Control -> Ambient en sd_sp
  levels(sd_sp[[cond_col]]) <- gsub("Control", "Ambient",
                                    levels(sd_sp[[cond_col]]))
  
  # Columna de nombre fijo para las fórmulas de adonis2. adonis2 falla con
  # "object 'sd_sp' not found" si la fórmula usa sd_sp[[loc_col]] directamente
  # y la llamada está dentro de una función — necesita un nombre de columna
  # literal, no una expresión con [[ ]].
  sd_sp$Locality_mod <- sd_sp[[loc_col]]
  
  # ---- PERMANOVA: Jaccard ~ Locality --------------------------
  set.seed(42)
  perm_res <- adonis2(
    jac ~ Locality_mod,
    data         = sd_sp,
    permutations = 999,
    by           = "term"
  )
  cat("\n--- PERMANOVA: Jaccard ~ Locality ---\n")
  print(perm_res)
  
  # ---- Pairwise PERMANOVA entre localidades (FDR) ----------------
  locs  <- levels(sd_sp[[loc_col]])
  pares <- combn(locs, 2, simplify = FALSE)
  
  pw_loc <- bind_rows(lapply(pares, function(par) {
    idx      <- sd_sp[[loc_col]] %in% par
    dist_sub <- as.dist(as.matrix(jac)[idx, idx])
    sd_sub   <- sd_sp[idx, ]
    set.seed(42)
    res <- adonis2(dist_sub ~ Locality_mod, data = sd_sub, permutations = 999)
    data.frame(
      Loc1 = par[1], Loc2 = par[2],
      R2   = round(res$R2[1], 3),
      F    = round(res$F[1],  3),
      p    = res$`Pr(>F)`[1]
    )
  }))
  pw_loc$p_adj <- round(p.adjust(pw_loc$p, method = "fdr"), 3)
  pw_loc$sig   <- cut(pw_loc$p_adj,
                      breaks = c(-Inf, .001, .01, .05, .1, Inf),
                      labels = c("***", "**", "*", ".", ""))
  
  # Annotate comparison type
  ctrl_locs <- c("La Bajita", "Charco Verde", "Cercana La Bajita")
  vent_locs <- c("El Faro", "La Bombilla")
  pw_loc$Comparison_type <- ifelse(
    (pw_loc$Loc1 %in% ctrl_locs & pw_loc$Loc2 %in% ctrl_locs),
    "Ambient-Ambient",
    ifelse(
      (pw_loc$Loc1 %in% vent_locs & pw_loc$Loc2 %in% vent_locs),
      "Vent-Vent",
      "Ambient-Vent"
    )
  )
  
  cat("\n--- Pairwise localities (FDR) ---\n")
  print(pw_loc[order(pw_loc$Comparison_type, pw_loc$p_adj), ])
  
  # ---- Betadisper ------------------------------------------------
  disp_cond     <- betadisper(jac, sd_sp[[cond_col]])
  disp_test     <- permutest(disp_cond, permutations = 999)
  cat("\n--- PERMDISP by Condition ---\n"); print(disp_test)
  
  disp_loc      <- betadisper(jac, sd_sp[[loc_col]])
  disp_loc_test <- permutest(disp_loc, permutations = 999)
  cat("\n--- PERMDISP by Locality ---\n"); print(disp_loc_test)
  
  # Guardar stats
  sink(file.path(OUTDIR, "Stats_beta_Jaccard.txt"), append = TRUE)
  cat("\n\n===", sp_name, "===\n")
  cat("--- PERMANOVA: Jaccard ~ Locality ---\n\n"); print(perm_res)
  cat("\n--- Pairwise localities (FDR) ---\n\n")
  print(pw_loc[order(pw_loc$Comparison_type, pw_loc$p_adj), ])
  cat("\n--- PERMDISP by Condition ---\n\n"); print(disp_test)
  cat("\n--- PERMDISP by Locality ---\n\n"); print(disp_loc_test)
  sink()
  
  write.csv(pw_loc,
            file.path(OUTDIR, paste0("Stats_beta_pairwise_", sp_label, ".csv")),
            row.names = FALSE)
  
  cat("Stats creadas:", sp_label, "\n")
  
  list(
    sp_label = sp_label, sp_name = sp_name,
    jac = jac, sd_sp = sd_sp, perm_res = perm_res,
    pw_loc = pw_loc, disp_test = disp_test, disp_loc_test = disp_loc_test
  )
}
 
 
## ---- 4. BETA DIVERSITY — PCoA (gráfico) ----
 
plot_pcoa_species <- function(res) {
  
  jac      <- res$jac
  sd_sp    <- res$sd_sp
  sp_name  <- res$sp_name
  perm_res <- res$perm_res
  
  pcoa_r  <- cmdscale(jac, k = 2, eig = TRUE)
  var_exp <- round(
    pcoa_r$eig[1:2] / sum(pcoa_r$eig[pcoa_r$eig > 0]) * 100, 1)
  
  pcoa_df <- as.data.frame(pcoa_r$points) %>%
    setNames(c("PCoA1", "PCoA2")) %>%
    rownames_to_column("Sample") %>%
    left_join(sd_sp %>% rownames_to_column("Sample"), by = "Sample") %>%
    mutate(
      Locality  = factor(.data[[loc_col]], levels = loc_order_all),
      Condition = factor(.data[[cond_col]], levels = c("Ambient", "Vent"))
    )
  
  p_val_loc <- perm_res$`Pr(>F)`[1]
  
  p_pcoa <- ggplot(pcoa_df,
                   aes(x = PCoA1, y = PCoA2,
                       colour = Locality, shape = Condition)) +
    geom_point(size = 3.5, alpha = 0.85) +
    stat_ellipse(aes(group = Condition, linetype = Condition),
                 type = "t", level = 0.95, linewidth = 0.8) +
    scale_colour_manual(values = pal_loc) +
    scale_shape_manual(values = c("Ambient" = 16, "Vent" = 17)) +
    scale_linetype_manual(values = c("Ambient" = "dashed", "Vent" = "solid")) +
    labs(title    = bquote(italic(.(sp_name))),
         subtitle = paste0("PERMANOVA p = ", p_val_loc,
                           " | PCoA1 = ", var_exp[1],
                           "% | PCoA2 = ", var_exp[2], "%"),
         x        = paste0("PCoA1 (", var_exp[1], "%)"),
         y        = paste0("PCoA2 (", var_exp[2], "%)"),
         colour   = "Locality",
         shape    = "Condition",
         linetype = "Condition") +
    theme_bw(base_size = 12) +
    theme(legend.position  = "right",
          plot.title       = element_text(face = "bold.italic", hjust = 0.5),
          plot.subtitle    = element_text(size = 9, colour = "grey40",
                                          hjust = 0.5))
  
  cat("Plot creado:", res$sp_label, "\n")
  
  invisible(p_pcoa)
}


## ---- 4. BETA DIVERSITY — Panel combinado (ambas especies) ----
 
build_pcoa_panel <- function(pcoa_plots) {
  
  panel_pcoa <- wrap_plots(pcoa_plots, ncol = 2) +
    plot_annotation(
      title   = "PCoA Jaccard (presence/absence) — KK samples | La Palma",
      tag_levels = "A"
    ) &
    theme(plot.tag = element_text(face = "bold", size = 13))
  
  print(panel_pcoa)
  ggsave(file.path(OUTDIR, "Fig_PCoA_panel.pdf"),
         plot = panel_pcoa, width = 16, height = 7)
  ggsave(file.path(OUTDIR, "Fig_PCoA_panel.png"),
         plot = panel_pcoa, width = 16, height = 7, dpi = 300)
  cat("Guardado: Fig_PCoA_panel\n")
  
  invisible(panel_pcoa)
}
 
res_arb2 <- run_beta_stats("Arbacia_lixula")
p_arb    <- plot_pcoa_species(res_arb2)
res_par2 <- run_beta_stats("Paracentrotus_lividus")
p_par    <- plot_pcoa_species(res_par2)
build_pcoa_panel(list(Arbacia_lixula = p_arb, Paracentrotus_lividus = p_par))
 
## ---- 4b. BETA DIVERSITY - PERMANOVA por Condición --------------
cat("\n=== 4b. Beta diversity — Jaccard PERMANOVA por Condición ===\n")
 
sink(file.path(OUTDIR, "Stats_beta_Jaccard_Condition.txt"))
cat("=== BETA DIVERSITY: Jaccard (presence/absence) ~ Condition ===\n")
cat("Model: adonis2(Jaccard ~ Condition)\n")
cat("Control renombrado a Ambient\n")
cat("Betadisper por Condition\n\n")
sink()
 
 
## ---- 4b. BETA DIVERSITY — PERMANOVA por Condición: ESTADÍSTICA ----
 
run_condition_stats <- function(sp_label) {
  
  ps_sp   <- ps_list[[sp_label]]
  sp_name <- sp_names[[sp_label]]
  cat("\n---", sp_name, "---\n")
  
  # Rarefy
  set.seed(42)
  ps_r <- rarefy_even_depth(ps_sp, sample.size = min_reads,
                            rngseed = 42, replace = FALSE,
                            trimOTUs = TRUE, verbose = FALSE)
  
  otu_r <- as(otu_table(ps_r), "matrix")
  if (!taxa_are_rows(ps_r)) otu_r <- t(otu_r)
  otu_t <- t(otu_r)
  
  jac <- vegdist(otu_t, method = "jaccard", binary = TRUE)
  
  sd_sp <- as.data.frame(as.matrix(sample_data(ps_r)))
  
  # Renombrar Control -> Ambient
  sd_sp[[cond_col]] <- gsub("Control", "Ambient", sd_sp[[cond_col]])
  sd_sp[[cond_col]] <- as.factor(sd_sp[[cond_col]])
  sd_sp[[loc_col]]  <- as.factor(sd_sp[[loc_col]])
  
  cat("Niveles Condition:", levels(sd_sp[[cond_col]]), "\n")
  cat("n Ambient:", sum(sd_sp[[cond_col]] == "Ambient"),
      "| n Vent:", sum(sd_sp[[cond_col]] == "Vent"), "\n")
  
  # Columnas de nombre fijo para las fórmulas — adonis2 falla con
  # "object 'sd_sp' not found" si la fórmula usa sd_sp[[col]] directamente
  # dentro de una función.
  sd_sp$Condition_mod <- sd_sp[[cond_col]]
  sd_sp$Locality_mod  <- sd_sp[[loc_col]]
  
  # ---- PERMANOVA: Jaccard ~ Condition -------------------------

  set.seed(42)
  perm_cond <- adonis2(
    jac ~ Condition_mod,
    data         = sd_sp,
    permutations = 999,
    by           = "term"
  )
  cat("\n--- PERMANOVA: Jaccard ~ Condition (DESCRIPTIVE — R2 only, no valid p) ---\n")
  cat("R2 =", round(perm_cond$R2[1], 3), "\n")
  
  # ---- PERMANOVA: Jaccard ~ Condition + Locality --------------

  set.seed(42)
  perm_cond_loc <- adonis2(
    jac ~ Condition_mod + Locality_mod,
    data         = sd_sp,
    permutations = 999,
    by           = "term"
  )
  cat("\n--- PERMANOVA: Jaccard ~ Condition + Locality (review before citing) ---\n")
  print(perm_cond_loc)
  
  # ---- Betadisper por Condición --------------------------------
  disp_cond     <- betadisper(jac, sd_sp[[cond_col]])
  disp_test     <- permutest(disp_cond, permutations = 999)
  cat("\n--- PERMDISP by Condition ---\n")
  print(disp_test)
  
  # ---- Guardar stats ------------------------------------------
  sink(file.path(OUTDIR, "Stats_beta_Jaccard_Condition.txt"), append = TRUE)
  cat("\n\n===", sp_name, "===\n")
  cat("n Ambient:", sum(sd_sp[[cond_col]] == "Ambient"),
      "| n Vent:", sum(sd_sp[[cond_col]] == "Vent"), "\n")
  cat("\n--- PERMANOVA: Jaccard ~ Condition (DESCRIPTIVE — R2 only) ---\n\n")
  cat("R2 =", round(perm_cond$R2[1], 3), "\n")
  cat("\n--- PERMANOVA: Jaccard ~ Condition + Locality (review before citing) ---\n\n")
  print(perm_cond_loc)
  cat("\n--- PERMDISP by Condition ---\n\n")
  print(disp_test)
  sink()
  
  cat("Stats creadas:", sp_label, "\n")
  
  list(
    sp_label = sp_label, sp_name = sp_name,
    jac = jac, sd_sp = sd_sp,
    perm_cond = perm_cond, perm_cond_loc = perm_cond_loc,
    disp_test = disp_test
  )
}
 
 
## ---- 4b. BETA DIVERSITY — PERMANOVA por Condición: PCoA (gráfico) ----
 
plot_pcoa_condition <- function(res) {
  
  jac       <- res$jac
  sd_sp     <- res$sd_sp
  sp_name   <- res$sp_name
  sp_label  <- res$sp_label
  perm_cond <- res$perm_cond
  
  pcoa_r  <- cmdscale(jac, k = 2, eig = TRUE)
  var_exp <- round(
    pcoa_r$eig[1:2] / sum(pcoa_r$eig[pcoa_r$eig > 0]) * 100, 1)
  
  pcoa_df <- as.data.frame(pcoa_r$points) %>%
    setNames(c("PCoA1", "PCoA2")) %>%
    rownames_to_column("Sample") %>%
    left_join(sd_sp %>% rownames_to_column("Sample"), by = "Sample") %>%
    mutate(
      Condition = factor(.data[[cond_col]], levels = c("Ambient", "Vent")),
      Locality  = factor(.data[[loc_col]],  levels = loc_order_all)
    )
  
  p_pcoa <- ggplot(pcoa_df,
                   aes(x = PCoA1, y = PCoA2,
                       colour = Condition, shape = Locality)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = Condition), level = 0.95,
                 linetype = "dashed", linewidth = 0.7) +
    scale_colour_manual(values = c("Ambient" = "#2166AC",
                                   "Vent"    = "#D73027")) +
    scale_shape_manual(values = c(16, 17, 15, 18, 8)) +
    labs(
      title    = paste("PCoA Jaccard (P/A) –", sp_name),
      subtitle = paste("Condition: R² =", round(perm_cond$R2[1], 3),
                       "(descriptive; no valid p — see comment above)"),
      x = paste0("PCoA1 (", var_exp[1], "%)"),
      y = paste0("PCoA2 (", var_exp[2], "%)"),
      colour = "Condition", shape = "Locality"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title    = element_text(face = "italic"),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      legend.position = "right"
    )
  
  ggsave(file.path(OUTDIR, paste0("Fig_PCoA_Condition_", sp_label, ".png")),
         plot = p_pcoa, width = 9, height = 7, dpi = 300)
  ggsave(file.path(OUTDIR, paste0("Fig_PCoA_Condition_", sp_label, ".pdf")),
         plot = p_pcoa, width = 9, height = 7)
  
  cat("Guardado: PCoA Condition –", sp_label, "\n")
  
  invisible(p_pcoa)
}
 
res_cond_arb <- run_condition_stats("Arbacia_lixula")
plot_pcoa_condition(res_cond_arb)
res_cond_par <- run_condition_stats("Paracentrotus_lividus")
plot_pcoa_condition(res_cond_par) 

## ---- 5. PERMANOVA POR PHYLUM Y FAMILY (FDR) -------------------
cat("\n=== 5. PERMANOVA por Phylum y Family ===\n")
 
## ---- Helper: safe PERMANOVA ------------------------------------
safe_permanova_taxon <- function(ps_sub, sp_name, taxrank, tax_i,
                                 loc_col, cond_col, nperm = 999) {
  
  if (nsamples(ps_sub) < 4) {
    cat("    [SKIP]", tax_i, ": n <4\n"); return(NULL)
  }
  
  sd_sub <- as.data.frame(sample_data(ps_sub))
  
  otu_mat <- as(otu_table(ps_sub), "matrix")
  if (!taxa_are_rows(ps_sub)) otu_mat <- t(otu_mat)
  otu_t <- t(otu_mat)
  
  # Filtrar filas (muestras) con suma 0
  keep_rows <- rowSums(otu_t) > 0
  if (sum(keep_rows) < 4) {
    cat("    [SKIP]", tax_i, ": <4 muestras con reads\n"); return(NULL)
  }
  otu_t  <- otu_t[keep_rows, , drop = FALSE]
  sd_sub <- sd_sub[keep_rows, , drop = FALSE]
  
  # Filtrar columnas (OTUs) con suma 0
  otu_t <- otu_t[, colSums(otu_t) > 0, drop = FALSE]
  if (ncol(otu_t) == 0) {
    cat("    [SKIP]", tax_i, ": sin OTUs\n"); return(NULL)
  }
  
  # Solo OTUs presentes en >= 2 muestras
  present_in <- colSums(otu_t > 0)
  otu_t <- otu_t[, present_in >= 2, drop = FALSE]
  if (ncol(otu_t) == 0) {
    cat("    [SKIP]", tax_i, ": no OTUs en >=2 muestras\n"); return(NULL)
  }
  
  # Mínimo 2 OTUs
  if (ncol(otu_t) < 2) {
    cat("    [SKIP]", tax_i, ": <2 OTUs\n"); return(NULL)
  }
  
  # Eliminar filas vacías tras filtrar OTUs
  keep_rows2 <- rowSums(otu_t) > 0
  if (sum(keep_rows2) < 4) {
    cat("    [SKIP]", tax_i, ": <4 muestras tras filtro OTUs\n"); return(NULL)
  }
  otu_t  <- otu_t[keep_rows2, , drop = FALSE]
  sd_sub <- sd_sub[keep_rows2, , drop = FALSE]
  
  # Jaccard P/A
  jac <- tryCatch(
    suppressWarnings(vegdist(otu_t, method = "jaccard", binary = TRUE)),
    error = function(e) { cat("    [ERR vegdist]", tax_i, "\n"); NULL }
  )
  if (is.null(jac) || any(is.nan(jac)) || all(jac == 0) || any(is.na(jac))) {
    cat("    [SKIP]", tax_i, ": jac degenerado\n"); return(NULL)
  }
  
  loc_v  <- as.factor(sd_sub[[loc_col]])
  cond_v <- as.factor(sd_sub[[cond_col]])
  levels(cond_v) <- gsub("Control", "Ambient", levels(cond_v))
  
  # PERMANOVA ~ Locality
  res <- tryCatch({
    set.seed(42)
    adonis2(jac ~ loc_v, permutations = nperm, by = "term")
  }, error = function(e) {
    cat("    [ERR adonis2]", tax_i, "\n"); NULL
  })
  if (is.null(res)) return(NULL)
  
  if (!is.na(res$F[1]) && abs(res$F[1]) > 1000) {
    cat("    [SKIP]", tax_i, ": F value aberrante (", res$F[1], ")\n")
    return(NULL)
  }
  
  # ---- Pairwise PERMANOVA entre localidades ----
  locs  <- levels(loc_v)
  pares <- combn(locs, 2, simplify = FALSE)
  
  pw <- bind_rows(lapply(pares, function(par) {
    idx <- loc_v %in% par
    if (sum(idx) < 4) return(NULL)
    dist_sub <- as.dist(as.matrix(jac)[idx, idx])
    loc_sub  <- droplevels(loc_v[idx])
    tryCatch({
      set.seed(42)
      r <- adonis2(dist_sub ~ loc_sub, permutations = nperm)
      data.frame(
        Loc1 = par[1], Loc2 = par[2],
        R2   = round(r$R2[1], 3),
        F    = round(r$F[1],  3),
        p    = r$`Pr(>F)`[1]
      )
    }, error = function(e) NULL)
  }))
  
  if (!is.null(pw) && nrow(pw) > 0) {
    pw$p_adj <- round(p.adjust(pw$p, method = "fdr"), 3)
    pw$sig   <- as.character(cut(pw$p_adj,
                                 breaks = c(-Inf, .001, .01, .05, .1, Inf),
                                 labels = c("***", "**", "*", ".", "")))
    ctrl_locs <- c("La Bajita", "Charco Verde", "Cercana La Bajita")
    vent_locs <- c("El Faro", "La Bombilla")
    pw$Comparison_type <- ifelse(
      pw$Loc1 %in% ctrl_locs & pw$Loc2 %in% ctrl_locs, "Ambient-Ambient",
      ifelse(pw$Loc1 %in% vent_locs & pw$Loc2 %in% vent_locs, "Vent-Vent",
             "Ambient-Vent"))
  } else {
    pw <- NULL
  }
  
  # ---- Betadisper ----
  disp_cond <- tryCatch({
    bd <- betadisper(jac, cond_v)
    pt <- permutest(bd, permutations = 999)
    list(F = round(pt$tab$F[1], 3), p = round(pt$tab$`Pr(>F)`[1], 3))
  }, error = function(e) list(F = NA, p = NA))
  
  disp_loc <- tryCatch({
    bd <- betadisper(jac, loc_v)
    pt <- permutest(bd, permutations = 999)
    list(F = round(pt$tab$F[1], 3), p = round(pt$tab$`Pr(>F)`[1], 3))
  }, error = function(e) list(F = NA, p = NA))
  
  cat("    [OK]", tax_i,
      "| n_OTUs:", ncol(otu_t),
      "| n_samp:", nrow(otu_t),
      "| p_Loc:", round(res$`Pr(>F)`[1], 3), "\n")
  
  list(
    summary = data.frame(
      Species          = sp_name,
      Rank             = taxrank,
      Taxon            = tax_i,
      n_OTUs           = ncol(otu_t),
      n_samples        = nrow(otu_t),
      R2_Loc           = round(res$R2[1], 3),
      F_Loc            = round(res$F[1],  3),
      p_Loc            = res$`Pr(>F)`[1],
      Betadisp_Cond_F  = disp_cond$F,
      Betadisp_Cond_p  = disp_cond$p,
      Betadisp_Loc_F   = disp_loc$F,
      Betadisp_Loc_p   = disp_loc$p,
      stringsAsFactors = FALSE
    ),
    pairwise = pw
  )
}
 

## ---- Helper: run para un rank (SIN rarefacción) ----------------
run_taxon_permanova <- function(ps_sp, sp_name, sp_label,
                                taxrank, outdir, nperm = 999) {
  
  cat("\n\n====", sp_name, "|", taxrank, "====\n")
  
  tax_df_full <- as.data.frame(tax_table(ps_sp))
  if (!(taxrank %in% colnames(tax_df_full))) {
    cat("  [SKIP] rank no encontrado:", taxrank, "\n"); return(NULL)
  }
  
  all_taxa <- na.omit(unique(as.character(tax_df_full[[taxrank]])))
  all_taxa <- all_taxa[all_taxa != ""]
  cat("  Taxa candidatos:", length(all_taxa), "\n")
  
  # SIN rarefacción — usar ps_sp directamente
  ps_work <- ps_sp
  
  res_list <- list()
  
  for (i in seq_along(all_taxa)) {
    tax_i <- all_taxa[i]
    cat("  [", i, "/", length(all_taxa), "]", tax_i, "\n")
    
    otu_ids <- rownames(tax_df_full)[
      !is.na(tax_df_full[[taxrank]]) & tax_df_full[[taxrank]] == tax_i
    ]
    if (length(otu_ids) == 0) next
    
    ps_sub <- tryCatch({
      ps_tmp <- prune_taxa(otu_ids, ps_work)
      ps_tmp <- prune_samples(sample_sums(ps_tmp) > 0, ps_tmp)
      if (ntaxa(ps_tmp) == 0 || nsamples(ps_tmp) == 0) NULL else ps_tmp
    }, error = function(e) NULL)
    
    if (is.null(ps_sub)) {
      cat("    [SKIP]", tax_i, ": objeto phyloseq vacío\n"); next
    }
    
    result <- safe_permanova_taxon(
      ps_sub   = ps_sub,
      sp_name  = sp_name,
      taxrank  = taxrank,
      tax_i    = tax_i,
      loc_col  = loc_col,
      cond_col = cond_col,
      nperm    = nperm
    )
    
    if (!is.null(result)) res_list[[tax_i]] <- result
    
    saveRDS(res_list,
            file.path(outdir,
                      paste0("TEMP_", taxrank, "_", sp_label, ".rds")))
  }
  
  if (length(res_list) == 0) {
    cat("  Sin resultados para", taxrank, "\n"); return(NULL)
  }
  
  # Tabla resumen con FDR
  res_df <- bind_rows(lapply(res_list, function(x) x$summary))
  res_df$p_adj <- round(p.adjust(res_df$p_Loc, method = "fdr"), 4)
  res_df$sig   <- as.character(cut(res_df$p_adj,
                                   breaks = c(-Inf, .001, .01, .05, .1, Inf),
                                   labels = c("***", "**", "*", ".", ""),
                                   right  = FALSE))
  
  cat("\n--- PERMANOVA por", taxrank, "(ordenado por p_Loc) ---\n")
  print(res_df[order(res_df$p_Loc), ])
  
  pw_all <- bind_rows(lapply(names(res_list), function(tx) {
    pw <- res_list[[tx]]$pairwise
    if (!is.null(pw)) pw$Taxon <- tx
    pw
  }))
  
  write.csv(res_df[order(res_df$p_Loc), ],
            file.path(outdir, paste0("Stats_PERMANOVA_", taxrank,
                                     "_", sp_label, ".csv")),
            row.names = FALSE)
  
  if (!is.null(pw_all) && nrow(pw_all) > 0) {
    write.csv(pw_all,
              file.path(outdir, paste0("Stats_pairwise_", taxrank,
                                       "_", sp_label, ".csv")),
              row.names = FALSE)
  }
  
  return(list(summary = res_df, pairwise = pw_all))
}
 
## ---- Ejecutar --------------------------------------------------
sink(file.path(OUTDIR, "Stats_taxon_PERMANOVA.txt"))
cat("=== PERMANOVA Jaccard P/A por Phylum y Family ===\n")
cat("Model: adonis2(Jaccard(P/A) ~ Locality)\n")
cat("Sin rarefacción — datos originales\n")
cat("Filtros: >=4 muestras (10%), >=2 OTUs presentes en >=2 muestras\n")
cat("FDR correction dentro de cada nivel taxonomico y especie\n")
cat("Betadisper por Condition y Locality\n")
cat("Fecha:", format(Sys.time()), "\n\n")
sink()
 
all_results_taxon <- list()
 
## ---- 5. PERMANOVA POR PHYLUM Y FAMILY — una especie a la vez ----
 
all_results_taxon <- list()
 
run_taxon_permanova_species <- function(sp_label) {
  
  ps_sp   <- ps_list[[sp_label]]
  sp_name <- sp_names[[sp_label]]
  
  cat("\n\n========================================\n")
  cat("Especie:", sp_name, "\n")
  cat("========================================\n")
  
  res_phy <- run_taxon_permanova(ps_sp, sp_name, sp_label,
                                 "Phylum", OUTDIR, nperm = 999)
  res_fam <- run_taxon_permanova(ps_sp, sp_name, sp_label,
                                 "Family", OUTDIR, nperm = 999)
  
  sink(file.path(OUTDIR, "Stats_taxon_PERMANOVA.txt"), append = TRUE)
  cat("\n\n===", sp_name, "===\n")
  cat("\n--- Phylum (ordenado por p_Loc) ---\n\n")
  if (!is.null(res_phy)) print(res_phy$summary[order(res_phy$summary$p_Loc), ])
  cat("\n--- Family (ordenado por p_Loc) ---\n\n")
  if (!is.null(res_fam)) print(res_fam$summary[order(res_fam$summary$p_Loc), ])
  sink()
  
  cat("Terminado:", sp_label, "\n")
  
  list(Phylum = res_phy, Family = res_fam)
}
 
 
## ---- 5. PERMANOVA POR PHYLUM Y FAMILY — combinar resultados ----
 
combine_taxon_permanova <- function() {
  
  all_tax_df <- bind_rows(lapply(all_results_taxon, function(x) {
    bind_rows(
      if (!is.null(x$Phylum)) x$Phylum$summary else NULL,
      if (!is.null(x$Family)) x$Family$summary else NULL
    )
  }))
  write.csv(all_tax_df,
            file.path(OUTDIR, "Stats_taxon_PERMANOVA_summary.csv"),
            row.names = FALSE)
  cat("\nGuardado: Stats_taxon_PERMANOVA_summary.csv\n")
  cat("\n=== SECCIÓN 5 COMPLETADA ===\n")
  
  invisible(all_tax_df)
}
 
 
all_results_taxon[["Arbacia_lixula"]]       <- run_taxon_permanova_species("Arbacia_lixula")
all_results_taxon[["Paracentrotus_lividus"]] <- run_taxon_permanova_species("Paracentrotus_lividus")
combine_taxon_permanova()

 
## ---- 7. INDVAL.g (PRESENCIA/AUSENCIA) -------------------------
cat("\n=== 7. IndVal.g (P/A) ===\n")
 
OUTDIR7 <- file.path(OUTDIR, "IndVal")
dir.create(OUTDIR7, showWarnings = FALSE)
 
## Paleta localidades
pal_loc_dot <- c(
  "Cercana La Bajita" = "#1B9E77",
  "Charco Verde"      = "#80B1D3",
  "El Faro"           = "#FDB462",
  "La Bajita"         = "#66C2A5",
  "La Bombilla"       = "#FB8072"
)
 
ctrl_locs <- c("La Bajita", "Charco Verde", "Cercana La Bajita")
vent_locs <- c("El Faro", "La Bombilla")
 
## ================================================================
## FUNCIÓN: correr IndVal y guardar tablas
## ================================================================
run_indval_PA <- function(ps_sp, sp_name, sp_label,
                          cond_col, loc_col, outdir, nperm = 9999) {
  
  cat("\n---", sp_name, "---\n")
  
  # Rarefy
  set.seed(42)
  ps_r <- rarefy_even_depth(ps_sp, sample.size = min_reads,
                            rngseed = 42, replace = FALSE,
                            trimOTUs = TRUE, verbose = FALSE)
  
  otu_r <- as(otu_table(ps_r), "matrix")
  if (!taxa_are_rows(ps_r)) otu_r <- t(otu_r)
  
  # P/A
  pa_mat <- (otu_r > 0) * 1
  n_pres <- rowSums(pa_mat)
  pa_mat <- pa_mat[n_pres >= 2, , drop = FALSE]
  cat("  OTUs en IndVal (presentes >=2 muestras):", nrow(pa_mat), "\n")
  
  # Metadata
  sd_sp <- as.data.frame(sample_data(ps_r))
  sd_sp[[cond_col]] <- gsub("Control", "Ambient", sd_sp[[cond_col]])
  tax_r  <- as.data.frame(tax_table(ps_r))
  loc_vec <- as.factor(sd_sp[[loc_col]])
  
  cat("  Localidades:", paste(levels(loc_vec), collapse = ", "), "\n")
  cat("  n Ambient:", sum(sd_sp[[cond_col]] == "Ambient"),
      "| n Vent:", sum(sd_sp[[cond_col]] == "Vent"), "\n")
  
  # multipatt IndVal.g P/A
  set.seed(42)
  ind_res <- multipatt(
    t(pa_mat),
    cluster = loc_vec,
    func    = "IndVal.g",
    control = how(nperm = nperm)
  )
  
  # Extraer A y B
  A_mat <- ind_res$A
  B_mat <- ind_res$B
  
  sign_df <- ind_res$sign %>%
    as.data.frame() %>%
    rownames_to_column("OTU")
  
  indval_AB <- bind_rows(lapply(seq_len(nrow(sign_df)), function(i) {
    otu_i <- sign_df$OTU[i]
    if (otu_i %in% rownames(A_mat)) {
      data.frame(OTU = otu_i,
                 A   = round(max(A_mat[otu_i, ], na.rm = TRUE), 3),
                 B   = round(max(B_mat[otu_i, ], na.rm = TRUE), 3))
    } else {
      data.frame(OTU = otu_i, A = NA, B = NA)
    }
  }))
  
  # Tabla resultados
  ind_df <- sign_df %>%
    mutate(
      p_adj = p.adjust(p.value, method = "fdr"),
      sig   = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        p_adj < 0.1   ~ "·",
        TRUE          ~ ""
      )
    ) %>%
    left_join(tax_r %>% rownames_to_column("OTU"), by = "OTU") %>%
    left_join(indval_AB, by = "OTU") %>%
    mutate(
      taxon_label = case_when(
        !is.na(Species) & Species != "" ~ Species,
        !is.na(Genus)   & Genus   != "" ~ paste0(Genus, " sp."),
        !is.na(Family)  & Family  != "" ~ paste0(Family, " (Family)"),
        !is.na(Phylum)  & Phylum  != "" ~ paste0(Phylum, " (Phylum)"),
        TRUE ~ OTU
      )
    )
  
  # Asignar assoc_group
  s_cols_present <- colnames(ind_df)[grepl("^s\\.", colnames(ind_df))]
  cat("  Columnas s.*:", paste(s_cols_present, collapse = ", "), "\n")
  
  s_ctrl <- paste0("s.", ctrl_locs)[paste0("s.", ctrl_locs) %in% s_cols_present]
  s_vent <- paste0("s.", vent_locs)[paste0("s.", vent_locs) %in% s_cols_present]
  
  cat("  s_ctrl:", paste(s_ctrl, collapse = ", "), "\n")
  cat("  s_vent:", paste(s_vent, collapse = ", "), "\n")
  
  ind_df <- ind_df %>%
    mutate(
      assoc_ctrl = rowSums(select(., all_of(s_ctrl)) > 0, na.rm = TRUE),
      assoc_vent = rowSums(select(., all_of(s_vent)) > 0, na.rm = TRUE),
      assoc_group = case_when(
        assoc_ctrl > 0 & assoc_vent > 0 ~ "Both",
        assoc_ctrl > 0                  ~ "Ambient",
        assoc_vent > 0                  ~ "Vent",
        TRUE                            ~ "Unknown"
      )
    )
  
  cat("  Distribución assoc_group (p<0.1):\n")
  print(table(ind_df$assoc_group[ind_df$p_adj < 0.1]))
  
  # Guardar tablas
  write.csv(ind_df[order(ind_df$p.value), ],
            file.path(outdir, paste0("Stats_IndVal_", sp_label, ".csv")),
            row.names = FALSE)
  
  write.csv(ind_df %>%
              filter(p_adj < 0.05,
                     assoc_group %in% c("Ambient", "Vent")) %>%
              mutate(IndVal = round(stat, 3)) %>%
              arrange(assoc_group, p_adj) %>%
              select(assoc_group, taxon_label, Phylum, Family, A, B, IndVal, p_adj),
            file.path(outdir, paste0("Stats_IndVal_sig05_", sp_label, ".csv")),
            row.names = FALSE)
  
  write.csv(ind_df %>%
              filter(p_adj < 0.10,
                     assoc_group %in% c("Ambient", "Vent")) %>%
              mutate(IndVal = round(stat, 3)) %>%
              arrange(assoc_group, p_adj) %>%
              select(assoc_group, taxon_label, Phylum, Family, A, B, IndVal, p_adj),
            file.path(outdir, paste0("Stats_IndVal_sig10_", sp_label, ".csv")),
            row.names = FALSE)
  
  cat("  Sig05 (Ambient+Vent):",
      sum(ind_df$p_adj < 0.05 & ind_df$assoc_group %in% c("Ambient","Vent"),
          na.rm = TRUE), "\n")
  cat("  Sig10 (Ambient+Vent):",
      sum(ind_df$p_adj < 0.10 & ind_df$assoc_group %in% c("Ambient","Vent"),
          na.rm = TRUE), "\n")
  cat("  Excluidos Both/Unknown (p<0.1):",
      sum(ind_df$p_adj < 0.10 & ind_df$assoc_group %in% c("Both","Unknown"),
          na.rm = TRUE), "\n")
  
  return(list(ind_df = ind_df, pa_mat = pa_mat, sd_sp = sd_sp, loc_vec = loc_vec))
}
 
## ================================================================
## FUNCIÓN: dotplot mejorado
## ================================================================
make_indval_dotplot <- function(res, ps_sp, sp_name, sp_label, outdir) {
  
  ind_df  <- res$ind_df
  pa_mat  <- res$pa_mat
  sd_sp   <- res$sd_sp
  loc_vec <- res$loc_vec
  locs_sp <- loc_order_all[loc_order_all %in% levels(loc_vec)]
  
  # Solo p_adj < 0.05, Ambient o Vent
  ind_sig <- ind_df %>%
    filter(p_adj < 0.05,
           assoc_group %in% c("Ambient", "Vent"))
  
  if (nrow(ind_sig) == 0) {
    cat("  Sin taxa sig05 para dotplot\n")
    return(NULL)
  }
  
  # Prevalencia por localidad
  prev_df <- lapply(locs_sp, function(loc) {
    samp_loc <- rownames(sd_sp)[sd_sp[[loc_col]] == loc]
    samp_loc <- samp_loc[samp_loc %in% colnames(pa_mat)]
    if (length(samp_loc) == 0) return(NULL)
    data.frame(
      OTU        = rownames(pa_mat),
      Locality   = loc,
      prevalence = rowMeans(pa_mat[, samp_loc, drop = FALSE]) * 100,
      n_samples  = length(samp_loc)
    )
  })
  
  prev_df <- bind_rows(prev_df) %>%
    filter(OTU %in% ind_sig$OTU) %>%
    left_join(ind_sig %>% select(OTU, taxon_label, Phylum,
                                 stat, p_adj, assoc_group),
              by = "OTU") %>%
    filter(!is.na(taxon_label)) %>%
    mutate(
      sig_plot = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE          ~ ""
      )
    )
  
  # Orden taxa
  phylum_ord <- prev_df %>%
    group_by(Phylum) %>%
    summarise(tot = sum(prevalence), .groups = "drop") %>%
    arrange(desc(tot)) %>%
    pull(Phylum)
  
  taxa_ord <- prev_df %>%
    distinct(taxon_label, Phylum, stat, assoc_group) %>%
    mutate(
      Phylum     = factor(Phylum, levels = phylum_ord),
      assoc_rank = ifelse(assoc_group == "Ambient", 1, 2)
    ) %>%
    arrange(assoc_rank, Phylum, desc(stat)) %>%
    pull(taxon_label) %>%
    unique()
  
  prev_df <- prev_df %>%
    mutate(
      taxon_label = factor(taxon_label, levels = taxa_ord),
      Phylum      = factor(Phylum, levels = phylum_ord),
      Locality    = factor(Locality, levels = locs_sp)
    )
  
  n_ctrl <- sum(locs_sp %in% ctrl_locs)
  n_locs <- length(locs_sp)
  n_taxa <- length(taxa_ord)
  pal_dyn <- pal_loc_dot[names(pal_loc_dot) %in% locs_sp]
  
  phy_present <- levels(prev_df$Phylum)
  pal_phy <- pal_phyla[names(pal_phyla) %in% phy_present]
  missing <- setdiff(phy_present, names(pal_phy))
  if (length(missing) > 0)
    pal_phy <- c(pal_phy,
                 setNames(colorRampPalette(c("#AAAAAA","#444444"))(length(missing)),
                          missing))
  
  loc_labels_named <- c(
    "La Bajita"         = "La Bajita\n(Ambient)",
    "Charco Verde"      = "Charco Verde\n(Ambient)",
    "Cercana La Bajita" = "Cercana La Bajita\n(Ambient)",
    "El Faro"           = "El Faro\n(Vent)",
    "La Bombilla"       = "La Bombilla\n(Vent)"
  )
  
  p_dot <- ggplot(prev_df,
                  aes(x = Locality, y = taxon_label)) +
    
    # Fondos ambient/vent
    annotate("rect", xmin = 0.5, xmax = n_ctrl + 0.5,
             ymin = -Inf, ymax = Inf,
             fill = "#DDEEFF", alpha = 0.15) +
    annotate("rect", xmin = n_ctrl + 0.5, xmax = n_locs + 0.5,
             ymin = -Inf, ymax = Inf,
             fill = "#FFDDDD", alpha = 0.15) +
    
    # Puntos
    geom_point(aes(size = prevalence, colour = Locality),
               alpha = 0.82, na.rm = TRUE) +
    
    # Asteriscos
    geom_text(
      data = prev_df %>%
        filter(sig_plot != "") %>%
        filter(
          (assoc_group == "Ambient" & Locality %in% ctrl_locs) |
            (assoc_group == "Vent"    & Locality %in% vent_locs)
        ) %>%
        group_by(taxon_label) %>%
        slice_max(prevalence, n = 1, with_ties = FALSE) %>%
        ungroup(),
      aes(label = sig_plot),
      nudge_x = 0.42, size = 3.2, colour = "grey20", fontface = "bold"
    ) +
    
    # Línea divisoria
    geom_vline(xintercept = n_ctrl + 0.5,
               linetype = "dashed", colour = "grey35", linewidth = 0.65) +
    
    facet_grid(Phylum ~ ., scales = "free_y", space = "free_y",
               switch = "y") +
    
    scale_size_continuous(
      range  = c(1, 13),
      limits = c(0, 100),
      breaks = c(20, 40, 60, 80, 100),
      name   = "Prevalence\n(% samples)"
    ) +
    scale_colour_manual(values = pal_dyn, guide = "none") +
    scale_x_discrete(
      labels = loc_labels_named[locs_sp],
      expand = expansion(add = c(0.7, 1.1))
    ) +
    scale_y_discrete(expand = expansion(add = c(0.6, 0.6))) +
    
    # Etiquetas Ambient/Vent UNA SOLA VEZ
    annotation_custom(
      grob = grid::textGrob(
        "Ambient",
        x = unit((n_ctrl / 2 + 0.5) / n_locs, "npc"),
        gp = grid::gpar(col = "#2166AC", fontface = "bold", fontsize = 10)
      ),
      ymin = Inf, ymax = Inf
    ) +
    annotation_custom(
      grob = grid::textGrob(
        "Vent",
        x = unit((n_ctrl + (n_locs - n_ctrl) / 2 + 0.5) / n_locs, "npc"),
        gp = grid::gpar(col = "#D73027", fontface = "bold", fontsize = 10)
      ),
      ymin = Inf, ymax = Inf
    ) +
    
    coord_cartesian(clip = "off") +
    
    labs(
      title    = bquote(italic(.(sp_name))),
      subtitle = "Indicator taxa (IndVal.g, P/A | FDR-adjusted, p < 0.05)\ndot size = prevalence (% samples) | * p<0.05 | ** p<0.01",
      x = NULL, y = NULL
    ) +
    
    theme_bw(base_size = 11) +
    theme(
      plot.title         = element_text(hjust = 0.5, face = "bold.italic",
                                        size = 13),
      plot.subtitle      = element_text(hjust = 0.5, size = 8.5,
                                        colour = "grey40"),
      axis.text.x        = element_text(size = 8.5, angle = 15,
                                        hjust = 1, vjust = 1),
      axis.text.y        = element_text(face = "italic", size = 8.5),
      strip.text.y.left  = element_text(angle = 0, face = "bold",
                                        size = 9, hjust = 1),
      strip.placement    = "outside",
      strip.background   = element_rect(fill = "#F0F0F0", colour = "grey70"),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      panel.spacing      = unit(0.3, "lines"),
      legend.position    = "right",
      legend.key.size    = unit(0.8, "lines"),
      plot.margin        = margin(t = 40, r = 10, b = 10, l = 10)
    )
  
  h_plot <- max(n_taxa * 0.45 + 3, 7)
  
  ggsave(file.path(outdir, paste0("Fig_IndVal_dotplot_", sp_label, ".pdf")),
         plot = p_dot, width = 10, height = h_plot)
  ggsave(file.path(outdir, paste0("Fig_IndVal_dotplot_", sp_label, ".png")),
         plot = p_dot, width = 10, height = h_plot, dpi = 300)
  
  cat("  Guardado dotplot –", sp_label, "\n")
  return(p_dot)
}
 
## ================================================================
## EJECUTAR
## ================================================================
res_arb <- run_indval_PA(
  ps_sp    = ps_arb,
  sp_name  = "Arbacia lixula",
  sp_label = "Arbacia_lixula",
  cond_col = cond_col,
  loc_col  = loc_col,
  outdir   = OUTDIR7,
  nperm    = 99990
)
 
res_par <- run_indval_PA(
  ps_sp    = ps_par,
  sp_name  = "Paracentrotus lividus",
  sp_label = "Paracentrotus_lividus",
  cond_col = cond_col,
  loc_col  = loc_col,
  outdir   = OUTDIR7,
  nperm    = 99990
)
 
# Guardar ind_df para uso posterior
ind_arb <- res_arb$ind_df
ind_par <- res_par$ind_df
 
# Dotplots
p_arb <- make_indval_dotplot(
  res      = res_arb,
  ps_sp    = ps_arb,
  sp_name  = "Arbacia lixula",
  sp_label = "Arbacia_lixula",
  outdir   = OUTDIR7
)
 
p_par <- make_indval_dotplot(
  res      = res_par,
  ps_sp    = ps_par,
  sp_name  = "Paracentrotus lividus",
  sp_label = "Paracentrotus_lividus",
  outdir   = OUTDIR7
)
 
cat("\n=== INDVAL P/A COMPLETADO ===\n")
cat("Outputs en:", OUTDIR7, "\n")
 
cat("\nArbacia lixula:\n")
cat("  Sig05:", sum(ind_arb$p_adj < 0.06 &
                      ind_arb$assoc_group %in% c("Ambient","Vent"), na.rm=TRUE), "\n")
cat("  Sig10:", sum(ind_arb$p_adj < 0.10 &
                      ind_arb$assoc_group %in% c("Ambient","Vent"), na.rm=TRUE), "\n")
 
cat("\nParacentrotus lividus:\n")
cat("  Sig05:", sum(ind_par$p_adj < 0.05 &
                      ind_par$assoc_group %in% c("Ambient","Vent"), na.rm=TRUE), "\n")
cat("  Sig10:", sum(ind_par$p_adj < 0.10 &
                      ind_par$assoc_group %in% c("Ambient","Vent"), na.rm=TRUE), "\n")
 

## ================================================================
