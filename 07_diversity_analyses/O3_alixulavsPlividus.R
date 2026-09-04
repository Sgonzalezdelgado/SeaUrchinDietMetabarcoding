### ══════════════════════════════════════════════════════════════════════════
### OBJECTIVE 3 — Interspecific dietary overlap, A. lixula vs P. lividus, 
### across localities
### ══════════════════════════════════════════════════════════════════════════

library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(RColorBrewer)
library(patchwork)

setwd(".")  # set to your working directory containing the final curated phyloseq object

OUTDIR <- "Objetivo3"
dir.create(OUTDIR, showWarnings = FALSE)

## ---- 1. CARGA Y FILTROS -----------------------------------------------------
cat("\n=== 1. Carga y filtros ===\n")

ps_full <- readRDS("phyloseq_OTU97_LULU_e20_EUK_KKE_curated_speciesGlom.rds")

sd_all <- as.data.frame(sample_data(ps_full))

# Detectar columnas automáticamente
kk_col   <- c("tissue","sample_type","matrix","Material","type")[
  c("tissue","sample_type","matrix","Material","type") %in% colnames(sd_all)][1]
cond_col <- c("condition","CONDICION","Condition","treatment")[
  c("condition","CONDICION","Condition","treatment") %in% colnames(sd_all)][1]
loc_col  <- c("locality","Localidad","location","Site")[
  c("locality","Localidad","location","Site") %in% colnames(sd_all)][1]
sp_col   <- c("host_species","ESPECIE","species","Erizo")[
  c("host_species","ESPECIE","species","Erizo") %in% colnames(sd_all)][1]

cat("Columnas — KK:", kk_col, "| Cond:", cond_col,
    "| Loc:", loc_col, "| Sp:", sp_col, "\n")

# ---- Selección de muestras --------------------------------------------------
# Solo KK + solo erizos (excluye agua W, roca R y sus NAs en sp_col)
ps <- subset_samples(ps_full, get(kk_col) == "KK")
ps <- subset_samples(ps,
                     get(sp_col) %in% c("Arbacia lixula", "Paracentrotus lividus"))

# Excluir Blanes, Los Porretos (fuera del diseño / n insuficiente)
ps <- subset_samples(ps,
                     !grepl("Blanes|Los Porretos|LaPorretos", get(loc_col),
                            ignore.case = TRUE))

# ---- Eliminar host (erizo) y Echinodermata de la tabla de OTUs --------------
tax_df   <- as.data.frame(tax_table(ps))
host_idx <- grepl("Arbacia|Paracentrotus", tax_df$Species, ignore.case = TRUE) |
  (!is.na(tax_df$Phylum) & tax_df$Phylum == "Echinodermata")

if (sum(host_idx) > 0) {
  ps <- prune_taxa(!host_idx, ps)
  cat("OTUs eliminados (host + Echinodermata):", sum(host_idx), "\n")
} else {
  cat("Sin OTUs host/Echinodermata en dataset.\n")
}

ps <- prune_taxa(taxa_sums(ps) > 0, ps)

# Renombrar Transition Faro -> Los Porretos  y  Transition -> Vent
# (por si quedó alguna etiqueta residual)
sd_fix <- as.data.frame(sample_data(ps))
sd_fix[[loc_col]]  <- gsub("Transition Faro|Transicion Faro",
                           "Los Porretos", sd_fix[[loc_col]])
sd_fix[[cond_col]] <- gsub("Transition|Transicion", "Vent", sd_fix[[cond_col]])
sample_data(ps) <- sample_data(sd_fix)

# Segunda pasada por si el rename creó Los Porretos nuevos
ps <- subset_samples(ps,
                     !grepl("Los Porretos|LaPorretos", get(loc_col), ignore.case = TRUE))
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps <- prune_samples(sample_sums(ps) > 0, ps)

cat("\n=== Objeto final tras filtros ===\n"); print(ps)

sd_ps <- as.data.frame(sample_data(ps))
cat("\nDistribucion Localidad x Especie:\n")
print(table(sd_ps[[loc_col]], sd_ps[[sp_col]]))

saveRDS(ps, file.path(OUTDIR, "ps_KK_LaPalma_bentonica.rds"))

## ---- PALETAS ----------------------------------------------------------------
sp_colors <- c(
  "Arbacia lixula"        = "#9B6FC4",
  "Paracentrotus lividus" = "#61AA65"
)

loc_levels_full <- c("La Bajita", "Cercana La Bajita", "Charco Verde",
                     "El Faro", "La Bombilla")
loc_ambient <- c("La Bajita", "Cercana La Bajita", "Charco Verde")
loc_vent    <- c("El Faro", "La Bombilla")

loc_palette <- c(
  "La Bajita"         = "#2166AC",
  "Cercana La Bajita" = "#4393C3",
  "Charco Verde"      = "#92C5DE",
  "El Faro"           = "#D73027",
  "La Bombilla"       = "#F46D43"
)

## ---- 2. RAREFACCIÓN ---------------------------------------------------------
cat("\n=== 2. Rarefaccion ===\n")

reads_per_sample <- sort(sample_sums(ps))
cat("Min:", min(reads_per_sample), "| Max:", max(reads_per_sample),
    "| Median:", median(reads_per_sample), "\n")
cat("Muestras < 5000 reads:", sum(reads_per_sample < 5000), "\n")

set.seed(42)
RAREDEPTH <- 5000
ps_rare <- rarefy_even_depth(ps, sample.size = RAREDEPTH,
                             rngseed = 42, replace = FALSE, trimOTUs = TRUE)
cat("Tras rarefaccion a", RAREDEPTH, "reads:\n"); print(ps_rare)

sd_rare <- as.data.frame(sample_data(ps_rare))

# Filtrar loc_levels a los que realmente existen en los datos rarefactados
loc_levels <- loc_levels_full[loc_levels_full %in% unique(sd_rare[[loc_col]])]
n_ambient  <- sum(loc_levels %in% loc_ambient)
loc_palette <- loc_palette[names(loc_palette) %in% loc_levels]

cat("Localidades (orden final):", paste(loc_levels, collapse = " | "), "\n")

## ---- 3. ALPHA DIVERSITY (Richness) — construir datos ------------------------

build_alpha_df <- function() {
  
  cat("\n=== 3. Alpha diversity ===\n")
  
  alpha_df <- estimate_richness(ps_rare, measures = "Observed")
  alpha_df$Sample <- rownames(alpha_df)
  
  sd_rare_df         <- as.data.frame(as.matrix(sample_data(ps_rare)))
  sd_rare_df$Sample  <- rownames(sd_rare_df)
  alpha_df           <- merge(alpha_df, sd_rare_df, by = "Sample", all.x = TRUE)
  alpha_df[[loc_col]] <- factor(alpha_df[[loc_col]], levels = loc_levels)
  
  # Eliminar NAs residuales (Blanes u otros excluidos)
  alpha_df <- alpha_df[!is.na(alpha_df[[loc_col]]) &
                         !is.na(alpha_df[[sp_col]]), ]
  
  write.csv(alpha_df[, c("Sample", "Observed", loc_col, sp_col)],
            file.path(OUTDIR, "alpha_diversity_table.csv"), row.names = FALSE)
  
  alpha_df
}


## ---- 3. ALPHA DIVERSITY — gráfico -------------------------------------------

plot_alpha_richness <- function(alpha_df) {
  
  p_rich <- ggplot(alpha_df,
                   aes(x = .data[[loc_col]], y = Observed,
                       fill = .data[[sp_col]])) +
    geom_boxplot(outlier.size = 1, alpha = 0.85,
                 position = position_dodge(0.8)) +
    geom_jitter(aes(color = .data[[sp_col]]),
                position = position_jitterdodge(jitter.width = 0.15,
                                                dodge.width  = 0.8),
                size = 1.5, alpha = 0.6) +
    geom_vline(xintercept = n_ambient + 0.5,
               linetype = "dashed", colour = "grey40", linewidth = 0.7) +
    annotate("text", x = n_ambient / 2 + 0.5, y = Inf,
             label = "Ambient", vjust = 1.5, colour = "#2166AC",
             fontface = "bold", size = 3.5) +
    annotate("text",
             x = n_ambient + (length(loc_levels) - n_ambient) / 2 + 0.5,
             y = Inf, label = "Vent", vjust = 1.5, colour = "#D73027",
             fontface = "bold", size = 3.5) +
    scale_fill_manual(values  = sp_colors) +
    scale_color_manual(values = sp_colors) +
    scale_x_discrete(drop = TRUE) +
    labs(title = "Alpha diversity (Richness) — KK bentonicas",
         x = NULL, y = "No. OTUs observados",
         fill = "Species", color = "Species") +
    theme_bw(base_size = 12) +
    theme(axis.text.x  = element_text(angle = 30, hjust = 1),
          legend.text   = element_text(face = "italic"),
          plot.title    = element_text(face = "bold"))
  
  ggsave(file.path(OUTDIR, "alpha_richness_locality.pdf"),
         p_rich, width = 9, height = 5)
  cat("Guardado: alpha_richness_locality.pdf\n")
  
  invisible(p_rich)
}


## ---- 3b. LM Richness por localidad ------------------------------------------

run_alpha_lm_by_locality <- function(alpha_df) {
  
  cat("\n--- LM Richness: A.lixula vs P.lividus por localidad ---\n")
  
  alpha_loc_res <- data.frame()
  
  for (loc in loc_levels) {
    sub  <- alpha_df[alpha_df[[loc_col]] == loc, ]
    n_sp <- table(sub[[sp_col]])
    if (length(n_sp) < 2 || any(n_sp < 2)) {
      cat("  Saltando", loc, "- pocas muestras\n"); next
    }
    
    # Columna de nombre fijo — lm() con get(sp_col) dentro de una función
    # puede fallar de la misma forma que adonis2; más seguro usar un nombre
    # literal en la fórmula.
    sub$Species_mod <- factor(sub[[sp_col]])
    lm_fit  <- lm(Observed ~ Species_mod, data = sub)
    lm_sum  <- summary(lm_fit)
    coef_tab <- lm_sum$coefficients
    
    beta_par <- if (nrow(coef_tab) >= 2) round(coef_tab[2, "Estimate"],   2) else NA
    se_par   <- if (nrow(coef_tab) >= 2) round(coef_tab[2, "Std. Error"], 2) else NA
    t_stat   <- if (nrow(coef_tab) >= 2) round(coef_tab[2, "t value"],    3) else NA
    p_raw    <- if (nrow(coef_tab) >= 2) coef_tab[2, "Pr(>|t|)"]             else NA
    
    arb <- sub$Observed[sub[[sp_col]] == "Arbacia lixula"]
    par <- sub$Observed[sub[[sp_col]] == "Paracentrotus lividus"]
    
    alpha_loc_res <- rbind(alpha_loc_res, data.frame(
      Locality        = loc,
      n_Arbacia       = length(arb),
      mean_Arbacia    = round(mean(arb), 1),
      sd_Arbacia      = round(sd(arb),   1),
      n_Paracentrotus = length(par),
      mean_Para       = round(mean(par), 1),
      sd_Para         = round(sd(par),   1),
      beta            = beta_par,
      SE              = se_par,
      t_stat          = t_stat,
      p_value         = round(p_raw, 4)
    ))
  }
  
  alpha_loc_res$p_adj <- round(p.adjust(alpha_loc_res$p_value, method = "BH"), 4)
  alpha_loc_res$sig   <- ifelse(alpha_loc_res$p_adj < 0.001, "***",
                                ifelse(alpha_loc_res$p_adj < 0.01,  "**",
                                       ifelse(alpha_loc_res$p_adj < 0.05,  "*",
                                              ifelse(alpha_loc_res$p_adj < 0.1,   ".", ""))))
  
  cat("\nRichness LM por especie y localidad:\n"); print(alpha_loc_res)
  write.csv(alpha_loc_res,
            file.path(OUTDIR, "LM_richness_especie_por_localidad.csv"),
            row.names = FALSE)
  
  alpha_loc_res
}

## ---- 4. BETA DIVERSITY - JACCARD P/A — construir matrices -------------------

build_beta_matrices <- function() {
  
  cat("\n=== 4. Beta diversity (Jaccard P/A) ===\n")
  
  otu_pa  <- as(otu_table(ps_rare), "matrix")
  if (!taxa_are_rows(ps_rare)) otu_pa <- t(otu_pa)
  otu_bin <- (otu_pa > 0) * 1
  jac     <- vegdist(t(otu_bin), method = "jaccard", binary = TRUE)
  
  perm_df <- as.data.frame(as.matrix(sample_data(ps_rare)))
  perm_df <- perm_df[rownames(as.matrix(jac)), , drop = FALSE]
  perm_df[[loc_col]] <- factor(perm_df[[loc_col]], levels = loc_levels)
  perm_df[[sp_col]]  <- factor(perm_df[[sp_col]])
  
  # Eliminar NAs (Blanes u otros no en loc_levels)
  cat("NAs en sp_col:", sum(is.na(perm_df[[sp_col]])), "\n")
  cat("NAs en loc_col:", sum(is.na(perm_df[[loc_col]])), "\n")
  
  keep_rows <- !is.na(perm_df[[sp_col]]) & !is.na(perm_df[[loc_col]])
  perm_df   <- perm_df[keep_rows, , drop = FALSE]
  jac       <- as.dist(as.matrix(jac)[keep_rows, keep_rows])
  
  cat("Muestras tras eliminar NAs:", nrow(perm_df), "\n")
  
  perm_df$Species_mod  <- perm_df[[sp_col]]
  perm_df$Locality_mod <- perm_df[[loc_col]]
  
  list(jac = jac, perm_df = perm_df)
}


## ---- 4a. PERMANOVA global ---------------------------------------------------

run_permanova_global <- function(beta) {
  
  cat("\n--- PERMANOVA Jaccard ~ Especie * Localidad ---\n")
  
  set.seed(42)
  perm_int <- adonis2(
    jac ~ Species_mod * Locality_mod,
    data = beta$perm_df, permutations = 999, by = "terms"
  )
  rownames(perm_int) <- gsub("Species_mod",  "Species",  rownames(perm_int))
  rownames(perm_int) <- gsub("Locality_mod", "Locality", rownames(perm_int))
  print(perm_int)
  capture.output(perm_int) |>
    writeLines(file.path(OUTDIR, "PERMANOVA_Species_x_Locality.txt"))
  
  perm_int
}


## ---- 4b. Beta-dispersion -----------------------------------------------------

run_betadisper_global <- function(beta) {
  
  cat("\n--- Beta dispersion ~ Localidad ---\n")
  bd      <- betadisper(beta$jac, group = beta$perm_df[[loc_col]])
  bd_test <- permutest(bd, permutations = 999)
  print(bd_test)
  
  cat("\n--- Beta dispersion ~ Especie ---\n")
  bd_sp      <- betadisper(beta$jac, group = beta$perm_df[[sp_col]])
  bd_sp_test <- permutest(bd_sp, permutations = 999)
  print(bd_sp_test)
  
  capture.output({
    cat("=== Betadisper ~ Locality ===\n"); print(bd_test)
    cat("\n=== Betadisper ~ Species ===\n"); print(bd_sp_test)
  }) |> writeLines(file.path(OUTDIR, "betadisper_Loc_Sp.txt"))
  
  list(bd_test = bd_test, bd_sp_test = bd_sp_test)
}


## ---- 4c. PERMANOVA por localidad --------------------------------------------

run_permanova_by_locality <- function(beta) {
  
  cat("\n--- PERMANOVA: A.lixula vs P.lividus por localidad ---\n")
  
  perm_loc_res <- data.frame()
  jac     <- beta$jac
  perm_df <- beta$perm_df
  
  for (loc in loc_levels) {
    idx      <- perm_df[[loc_col]] == loc
    sub_df   <- perm_df[idx, , drop = FALSE]
    sub_dist <- as.dist(as.matrix(jac)[idx, idx])
    n_sp     <- table(sub_df[[sp_col]])
    
    if (length(n_sp) < 2 || any(n_sp < 2)) {
      cat("  Saltando", loc, "- pocas muestras\n"); next
    }
    
    set.seed(42)
    res    <- adonis2(sub_dist ~ Species_mod, data = sub_df, permutations = 999)
    bd_loc <- betadisper(sub_dist, group = sub_df[[sp_col]])
    bd_p   <- permutest(bd_loc, permutations = 999)$tab$`Pr(>F)`[1]
    
    perm_loc_res <- rbind(perm_loc_res, data.frame(
      Locality        = loc,
      n_Arbacia       = as.integer(n_sp["Arbacia lixula"]),
      n_Paracentrotus = as.integer(n_sp["Paracentrotus lividus"]),
      R2              = round(res$R2[1], 4),
      F_stat          = round(res$F[1],  4),
      p_value         = round(res$`Pr(>F)`[1], 4),
      betadisper_F    = round(permutest(bd_loc, permutations = 999)$tab$F[1], 3),
      betadisper_p    = round(bd_p, 4)
    ))
  }
  
  perm_loc_res$sig <- ifelse(perm_loc_res$p_value < 0.001, "***",
                             ifelse(perm_loc_res$p_value < 0.01,  "**",
                                    ifelse(perm_loc_res$p_value < 0.05,  "*",
                                           ifelse(perm_loc_res$p_value < 0.1,   ".", ""))))
  
  cat("\nResultados PERMANOVA por localidad:\n"); print(perm_loc_res)
  write.csv(perm_loc_res,
            file.path(OUTDIR, "PERMANOVA_por_localidad_especie.csv"),
            row.names = FALSE)
  
  perm_loc_res
}


## ---- 4d. TABLA UNIFICADA ----------------------------------------------------

build_unified_table <- function(perm_int, bd, perm_loc_res) {
  
  cat("\n--- Tabla unificada manuscrito ---\n")
  
  bd_test    <- bd$bd_test
  bd_sp_test <- bd$bd_sp_test
  
  global_tab <- data.frame(
    Section    = "PERMANOVA Global",
    Comparison = rownames(perm_int)[1:(nrow(perm_int)-2)],
    Df         = perm_int$Df[1:(nrow(perm_int)-2)],
    R2         = round(perm_int$R2[1:(nrow(perm_int)-2)], 3),
    F          = round(perm_int$F[1:(nrow(perm_int)-2)],  2),
    p          = perm_int$`Pr(>F)`[1:(nrow(perm_int)-2)]
  )
  
  betadisp_tab <- data.frame(
    Section    = "Betadisper Global",
    Comparison = c("~ Locality", "~ Species"),
    Df         = c(bd_test$tab$Df[1],    bd_sp_test$tab$Df[1]),
    R2         = NA,
    F          = round(c(bd_test$tab$F[1], bd_sp_test$tab$F[1]), 2),
    p          = c(bd_test$tab$`Pr(>F)`[1], bd_sp_test$tab$`Pr(>F)`[1])
  )
  
  pairwise_tab <- data.frame(
    Section    = "PERMANOVA por Localidad",
    Comparison = paste0("A. lixula vs P. lividus | ", perm_loc_res$Locality,
                        " (n=", perm_loc_res$n_Arbacia, " vs ",
                        perm_loc_res$n_Paracentrotus, ")"),
    Df         = 1,
    R2         = perm_loc_res$R2,
    F          = perm_loc_res$F_stat,
    p          = perm_loc_res$p_value
  )
  
  unified_tab <- rbind(global_tab, betadisp_tab, pairwise_tab)
  write.csv(unified_tab,
            file.path(OUTDIR, "PERMANOVA_tabla_unificada.csv"),
            row.names = FALSE)
  cat("Guardado: PERMANOVA_tabla_unificada.csv\n"); print(unified_tab)
  
  unified_tab
}


## ---- 4e. PCoA GLOBAL --------------------------------------------------------

plot_pcoa_global <- function(beta) {
  
  cat("\n--- PCoA global ---\n")
  
  jac     <- beta$jac
  perm_df <- beta$perm_df
  
  pcoa_res <- cmdscale(jac, k = 2, eig = TRUE)
  eig_vals <- pcoa_res$eig
  pct_var  <- round(100 * eig_vals / sum(eig_vals[eig_vals > 0]), 1)
  
  pcoa_df <- as.data.frame(pcoa_res$points) |>
    setNames(c("PCoA1", "PCoA2")) |>
    rownames_to_column("Sample") |>
    left_join(perm_df |> rownames_to_column("Sample"), by = "Sample")
  pcoa_df[[loc_col]] <- factor(pcoa_df[[loc_col]], levels = loc_levels)
  
  p_pcoa_global <- ggplot(pcoa_df,
                          aes(x = PCoA1, y = PCoA2,
                              color = .data[[loc_col]],
                              shape = .data[[sp_col]])) +
    geom_point(size = 3.5, alpha = 0.85) +
    stat_ellipse(aes(group = .data[[loc_col]]),
                 level = 0.75, linetype = "dashed", linewidth = 0.6) +
    scale_color_manual(values = loc_palette) +
    scale_shape_manual(values = c(16, 17),
                       labels = c("A. lixula", "P. lividus")) +
    labs(title    = "PCoA Jaccard — all localities",
         subtitle = "Fecal samples | benthic diet OTUs",
         x = paste0("PCoA1 (", pct_var[1], "%)"),
         y = paste0("PCoA2 (", pct_var[2], "%)"),
         color = "Locality", shape = "Species") +
    theme_bw(base_size = 12) +
    theme(legend.text = element_text(face = "italic"))
  
  ggsave(file.path(OUTDIR, "PCoA_global_locality_species.pdf"),
         p_pcoa_global, width = 9, height = 6)
  cat("Guardado: PCoA_global_locality_species.pdf\n")
  
  invisible(p_pcoa_global)
}


## ---- 4f. PCoA SEPARADOS por localidad — una localidad a la vez -------------

run_pcoa_one_locality <- function(beta, perm_loc_res, loc) {
  
  jac     <- beta$jac
  perm_df <- beta$perm_df
  sp_shapes <- c("Arbacia lixula" = 16, "Paracentrotus lividus" = 17)
  
  samps_loc <- rownames(perm_df)[as.character(perm_df[[loc_col]]) == loc]
  
  if (length(samps_loc) < 4) {
    cat("  Saltando", loc, "— solo", length(samps_loc), "muestras\n")
    return(NULL)
  }
  
  jac_loc  <- as.dist(as.matrix(jac)[samps_loc, samps_loc])
  pcoa_loc <- cmdscale(jac_loc, k = 2, eig = TRUE)
  eig_loc  <- pcoa_loc$eig
  pct_loc  <- round(100 * eig_loc / sum(eig_loc[eig_loc > 0]), 1)
  
  df_loc <- as.data.frame(pcoa_loc$points) |>
    setNames(c("PCoA1", "PCoA2")) |>
    rownames_to_column("Sample") |>
    left_join(perm_df |> rownames_to_column("Sample"), by = "Sample")
  df_loc[[sp_col]] <- factor(df_loc[[sp_col]])
  
  n_per_sp <- df_loc |> count(.data[[sp_col]]) |>
    filter(n >= 3) |> pull(.data[[sp_col]])
  
  perm_row <- perm_loc_res[perm_loc_res$Locality == loc, ]
  perm_sub <- if (nrow(perm_row) > 0) {
    paste0("R² = ", perm_row$R2, ", p = ", perm_row$p_value,
           " ", perm_row$sig)
  } else { "" }
  
  p_loc <- ggplot(df_loc,
                  aes(x = PCoA1, y = PCoA2,
                      color = .data[[sp_col]],
                      shape = .data[[sp_col]])) +
    geom_point(size = 4, alpha = 0.9) +
    {if (length(n_per_sp) >= 1)
      stat_ellipse(
        data      = df_loc |> filter(.data[[sp_col]] %in% n_per_sp),
        aes(group = .data[[sp_col]]),
        level = 0.75, linetype = "dashed", linewidth = 0.7)
    } +
    scale_color_manual(values = sp_colors, drop = FALSE,
                       labels = c("Arbacia lixula"        = "A. lixula",
                                  "Paracentrotus lividus" = "P. lividus")) +
    scale_shape_manual(values = sp_shapes, drop = FALSE,
                       labels = c("Arbacia lixula"        = "A. lixula",
                                  "Paracentrotus lividus" = "P. lividus")) +
    labs(title    = paste0(loc,
                           ifelse(loc %in% loc_vent, "  [Vent]", "  [Ambient]")),
         subtitle = paste0("n = ", length(samps_loc), " | ", perm_sub),
         x = paste0("PCoA1 (", pct_loc[1], "%)"),
         y = paste0("PCoA2 (", pct_loc[2], "%)"),
         color = NULL, shape = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          legend.position = "bottom",
          legend.text     = element_text(face = "italic"))
  
  safe_name <- gsub("[^A-Za-z0-9]", "_", loc)
  ggsave(file.path(OUTDIR, paste0("PCoA_", safe_name, ".pdf")),
         p_loc, width = 5.5, height = 5)
  cat(" ", loc, "— PCoA1:", pct_loc[1], "% | PCoA2:", pct_loc[2], "%\n")
  
  p_loc
}


## ---- 4f. Panel combinado (todas las localidades) ----------------------------

build_pcoa_locality_panel <- function(plot_list_pcoa) {
  
  plot_list_pcoa <- Filter(Negate(is.null), plot_list_pcoa)
  
  if (length(plot_list_pcoa) >= 2) {
    panel_pcoa <- wrap_plots(plot_list_pcoa, ncol = 3, guides = "collect") &
      theme(legend.position = "bottom")
    ggsave(file.path(OUTDIR, "PCoA_panel_localidades.pdf"),
           panel_pcoa, width = 16, height = 11)
    cat("Guardado: PCoA_panel_localidades.pdf\n")
    return(invisible(panel_pcoa))
  }
  
  cat("Menos de 2 localidades con panel válido — no se genera el panel.\n")
  invisible(NULL)
}

## ---- 5. SOLAPAMIENTO — helper por rango taxonómico (ya era función, sin cambios) ----

shared_by_rank <- function(ps_obj, rank, loc_levels, loc_col, sp_col) {
  
  ps_r      <- tax_glom(ps_obj, taxrank = rank, NArm = TRUE)
  otu_r     <- as(otu_table(ps_r), "matrix")
  if (!taxa_are_rows(ps_r)) otu_r <- t(otu_r)
  otu_r_bin <- (otu_r > 0) * 1
  tax_r     <- as.data.frame(tax_table(ps_r))
  sd_r      <- as.data.frame(as.matrix(sample_data(ps_obj)))
  res       <- data.frame()
  
  for (loc in loc_levels) {
    samps_arb <- intersect(
      rownames(sd_r)[sd_r[[loc_col]] == loc &
                       !is.na(sd_r[[sp_col]]) &
                       sd_r[[sp_col]] == "Arbacia lixula"],
      colnames(otu_r_bin))
    samps_par <- intersect(
      rownames(sd_r)[sd_r[[loc_col]] == loc &
                       !is.na(sd_r[[sp_col]]) &
                       sd_r[[sp_col]] == "Paracentrotus lividus"],
      colnames(otu_r_bin))
    
    if (length(samps_arb) == 0 | length(samps_par) == 0) next
    
    taxa_arb  <- rownames(otu_r_bin)[rowSums(otu_r_bin[, samps_arb, drop=FALSE]) > 0]
    taxa_par  <- rownames(otu_r_bin)[rowSums(otu_r_bin[, samps_par, drop=FALSE]) > 0]
    names_arb <- unique(tax_r[taxa_arb, rank])
    names_par <- unique(tax_r[taxa_par, rank])
    names_arb <- names_arb[!is.na(names_arb) & names_arb != ""]
    names_par <- names_par[!is.na(names_par) & names_par != ""]
    
    n_shared   <- length(intersect(names_arb, names_par))
    n_union    <- length(union(names_arb, names_par))
    n_excl_arb <- length(setdiff(names_arb, names_par))
    n_excl_par <- length(setdiff(names_par, names_arb))
    if (n_union == 0) next
    
    res <- rbind(res, data.frame(
      Localidad      = loc,
      Rank           = rank,
      N_excl_Arb    = n_excl_arb,
      N_compartidos  = n_shared,
      N_excl_Par    = n_excl_par,
      Total_union    = n_union,
      Pct_excl_Arb   = round(100 * n_excl_arb / n_union, 1),
      Pct_compartido = round(100 * n_shared    / n_union, 1),
      Pct_excl_Par   = round(100 * n_excl_par  / n_union, 1),
      Jaccard        = round(n_shared / n_union, 3)
    ))
  }
  return(res)
}


## ---- 5. SOLAPAMIENTO — correr los 4 rangos y guardar ------------------------

run_overlap_analysis <- function() {
  
  cat("\n=== 5. Solapamiento por rango taxonomico ===\n")
  
  rank_order <- c("Phylum", "Order", "Family", "Species")
  
  all_ranks_res <- do.call(rbind, lapply(rank_order, function(r) {
    cat("  Calculando rango:", r, "\n")
    shared_by_rank(ps_rare, r, loc_levels, loc_col, sp_col)
  }))
  
  write.csv(all_ranks_res,
            file.path(OUTDIR, "Solapamiento_por_rango_taxonomico.csv"),
            row.names = FALSE)
  
  all_ranks_res
}


## ---- 5. SOLAPAMIENTO — panel de gráficos por rango --------------------------

plot_overlap_panel <- function(all_ranks_res) {
  
  rank_order <- c("Phylum", "Order", "Family", "Species")
  
  pal_shared <- c(
    "Exclusivos A. lixula"  = "#9B6FC4",
    "Compartidos"           = "#AAAAAA",
    "Exclusivos P. lividus" = "#61AA65"
  )
  
  plot_list_ranks <- list()
  
  for (r in rank_order) {
    df_r <- all_ranks_res[all_ranks_res$Rank == r, ] |>
      select(Localidad, Pct_excl_Arb, Pct_compartido, Pct_excl_Par) |>
      pivot_longer(-Localidad, names_to = "Tipo", values_to = "Pct") |>
      mutate(
        Tipo = factor(Tipo,
                      levels = c("Pct_excl_Arb", "Pct_compartido", "Pct_excl_Par"),
                      labels = c("Exclusivos A. lixula", "Compartidos",
                                 "Exclusivos P. lividus")),
        Localidad = factor(Localidad, levels = loc_levels)
      ) |>
      filter(!is.na(Pct))
    
    if (nrow(df_r) == 0) next
    
    locs_con_datos <- unique(as.character(df_r$Localidad))
    n_ambient_r    <- sum(loc_ambient %in% locs_con_datos)
    
    p_r <- ggplot(df_r, aes(x = Localidad, y = Pct, fill = Tipo)) +
      geom_col(width = 0.7, color = "white", linewidth = 0.2, na.rm = TRUE) +
      geom_text(aes(label = ifelse(Pct > 5, paste0(Pct, "%"), "")),
                position = position_stack(vjust = 0.5),
                size = 2.8, color = "white", fontface = "bold", na.rm = TRUE) +
      annotate("segment",
               x = n_ambient_r + 0.5, xend = n_ambient_r + 0.5,
               y = 0, yend = 101,
               color = "grey30", linewidth = 0.7, linetype = "dashed") +
      scale_fill_manual(values = pal_shared) +
      scale_y_continuous(labels = function(x) paste0(x, "%"),
                         limits = c(0, 101), expand = c(0, 0)) +
      scale_x_discrete(drop = TRUE) +
      labs(title = r, x = NULL, y = "% taxa", fill = NULL) +
      theme_bw(base_size = 10) +
      theme(axis.text.x     = element_text(angle = 35, hjust = 1),
            legend.position = "none",
            plot.title      = element_text(face = "bold", size = 11))
    
    plot_list_ranks[[r]] <- p_r
  }
  
  panel_ranks <- wrap_plots(plot_list_ranks, ncol = 2, guides = "collect") +
    plot_annotation(
      title    = "Dietary overlap by taxonomic rank",
      subtitle = "% exclusive vs shared taxa | Ambient (left) → Vent (right)",
      theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                       plot.subtitle = element_text(size = 10, color = "grey40"))
    ) &
    theme(legend.position = "bottom")
  
  ggsave(file.path(OUTDIR, "Solapamiento_por_rango_taxonomico.pdf"),
         panel_ranks, width = 12, height = 14)
  ggsave(file.path(OUTDIR, "Solapamiento_por_rango_taxonomico.png"),
         panel_ranks, width = 12, height = 14, dpi = 300)
  cat("Guardado: Solapamiento_por_rango_taxonomico.pdf\n")
  
  invisible(panel_ranks)
}


## ---- 6. COMPOSICION TAXONOMICA POR LOCALIDAD (Phylum bar plot) --------------

plot_phylum_composition <- function() {
  
  cat("\n=== 6. Composicion Phylum x Localidad x Especie ===\n")
  
  ps_phylum <- tax_glom(ps_rare, taxrank = "Phylum", NArm = TRUE)
  
  df_phylum <- psmelt(ps_phylum) |>
    rename(Loc = all_of(loc_col), Sp = all_of(sp_col)) |>
    mutate(
      Phylum = ifelse(is.na(Phylum) | Phylum == "", "Unassigned", Phylum)
    ) |>
    # FIX CLAVE: filtrar Blanes y cualquier muestra cuya localidad no esté
    # en loc_levels (se convierten en NA al factorizar -> NA en el plot)
    filter(Loc %in% loc_levels,
           Sp  %in% c("Arbacia lixula", "Paracentrotus lividus")) |>
    mutate(Loc = factor(Loc, levels = loc_levels)) |>
    group_by(Sample, Loc, Sp) |>
    mutate(RelAbund = Abundance / sum(Abundance)) |>
    ungroup() |>
    filter(!is.nan(RelAbund))
  
  cat("Filas en df_phylum:", nrow(df_phylum), "\n")
  cat("NAs en Sp:", sum(is.na(df_phylum$Sp)), "\n")
  cat("NAs en Loc:", sum(is.na(df_phylum$Loc)), "\n")
  cat("Localidades presentes:", paste(unique(df_phylum$Loc), collapse = ", "), "\n")
  
  phylum_colors <- c(
    "Rhodophyta"        = "#E15759",
    "Ochrophyta"        = "#F28E2B",
    "Annelida"          = "#4E79A7",
    "Arthropoda"        = "#76B7B2",
    "Cnidaria"          = "#59A14F",
    "Mollusca"          = "#EDC948",
    "Porifera"          = "#B07AA1",
    "Platyhelminthes"   = "#FF9DA7",
    "Bacillariophyta"   = "#9C755F",
    "Chordata"          = "#BAB0AC",
    "Bryozoa"           = "#D4A6C8",
    "Chlorophyta"       = "#86BCB6",
    "Haptophyta"        = "#A0CBE8",
    "Nemertea"          = "#FFBE7D",
    "Nematoda"          = "#8CD17D",
    "Prasinodermophyta" = "#D4A6C8",
    "Rotifera"          = "#86BCB6",
    "Tardigrada"        = "#A0CBE8",
    "Unassigned"        = "#CCCCCC",
    "Other"             = "#C7C7C7"
  )
  
  # ---- Umbral por abundancia relativa media, calculado por Localidad x Especie ----
  THRESHOLD <- 0.01   # 1% en AL MENOS UNA localidad, por especie
  
  indval_sig_phyla <- c()   # p.ej. c("Annelida", "Arthropoda", "Rhodophyta")
  
  phyla_by_threshold <- df_phylum |>
    group_by(Sp, Loc, Phylum) |>
    summarise(MeanRel = mean(RelAbund, na.rm = TRUE), .groups = "drop") |>
    group_by(Sp, Phylum) |>
    summarise(MaxRel = max(MeanRel, na.rm = TRUE), .groups = "drop") |>
    filter(MaxRel > THRESHOLD) |>
    pull(Phylum) |>
    unique()
  
  top_phyla <- union(phyla_by_threshold, indval_sig_phyla)
  
  cat("\n--- Phyla mostrados individualmente (>", THRESHOLD * 100,
      "% en al menos una localidad, por especie) ---\n", sep = "")
  cat(paste(sort(top_phyla), collapse = ", "), "\n")
  
  todos_phyla <- unique(df_phylum$Phylum)
  cat("\n--- Phyla agrupados en 'Other' ---\n")
  cat(paste(sort(setdiff(todos_phyla, top_phyla)), collapse = ", "), "\n")
  
  df_sum_ph <- df_phylum |>
    mutate(Phylum_plot = ifelse(Phylum %in% top_phyla, Phylum, "Other")) |>
    group_by(Loc, Sp, Phylum_plot) |>
    summarise(MeanRel = mean(RelAbund, na.rm = TRUE), .groups = "drop") |>
    group_by(Loc, Sp) |>
    mutate(MeanRel = MeanRel / sum(MeanRel)) |>
    ungroup()
  
  # ---- Orden de apilado: algas (abajo) -> invertebrados (arriba) ----
  algae_order <- c("Bacillariophyta", "Ochrophyta", "Rhodophyta",
                   "Chlorophyta", "Haptophyta", "Prasinodermophyta")
  
  invert_order <- c("Annelida", "Arthropoda", "Mollusca", "Cnidaria",
                    "Porifera", "Platyhelminthes", "Chordata", "Bryozoa",
                    "Nemertea", "Nematoda", "Rotifera", "Tardigrada")
  
  present_levels <- unique(as.character(df_sum_ph$Phylum_plot))
  
  ordered_levels <- c(
    intersect(algae_order, present_levels),
    intersect(invert_order, present_levels),
    intersect(c("Unassigned", "Other"), present_levels)
  )
  
  faltantes <- setdiff(present_levels, ordered_levels)
  if (length(faltantes) > 0) {
    cat("AVISO: phyla no clasificados en algae_order/invert_order, añadidos al final:",
        paste(faltantes, collapse = ", "), "\n")
    ordered_levels <- c(ordered_levels, faltantes)
  }
  
  df_sum_ph <- df_sum_ph |>
    mutate(Phylum_plot = factor(Phylum_plot, levels = ordered_levels))
  
  locs_en_plot   <- levels(droplevels(df_sum_ph$Loc))
  n_ambient_plot <- sum(locs_en_plot %in% loc_ambient)
  
  p_phylum <- ggplot(df_sum_ph, aes(x = Loc, y = MeanRel, fill = Phylum_plot)) +
    geom_col(width = 0.8, color = "white", linewidth = 0.2,
             position = position_stack(reverse = TRUE)) +
    geom_vline(xintercept = n_ambient_plot + 0.5,
               linetype = "dashed", colour = "grey40", linewidth = 0.7) +
    facet_wrap(~ Sp, labeller = labeller(
      .default = function(x) gsub("Paracentrotus lividus", "P. lividus",
                                  gsub("Arbacia lixula", "A. lixula", x)))) +
    scale_fill_manual(values = phylum_colors, na.value = "#C7C7C7") +
    scale_y_continuous(labels = percent) +
    scale_x_discrete(drop = TRUE) +
    labs(x = NULL, y = "Mean relative abundance", fill = "Phylum") +
    theme_bw(base_size = 11) +
    theme(axis.text.x    = element_text(angle = 35, hjust = 1),
          strip.text      = element_text(face = "italic"),
          legend.position = "right")
  
  ggsave(file.path(OUTDIR, "Phylum_by_locality_species.pdf"),
         p_phylum, width = 13, height = 6)
  cat("Guardado: Phylum_by_locality_species.pdf\n")
  
  tabla_completa <- df_phylum |>
    group_by(Sp, Loc, Phylum) |>
    summarise(MeanRel = mean(RelAbund, na.rm = TRUE), .groups = "drop") |>
    group_by(Sp, Loc) |>
    mutate(Pct = round(MeanRel / sum(MeanRel) * 100, 2)) |>
    ungroup() |>
    select(Sp, Loc, Phylum, Pct) |>
    pivot_wider(names_from = Phylum, values_from = Pct, values_fill = 0) |>
    arrange(Sp, Loc)
  
  write.csv2(tabla_completa, file.path(OUTDIR, "Table_phylum_pct_completo.csv"),
             row.names = FALSE)
  cat("Guardado: Table_phylum_pct_completo.csv\n")
  
  list(plot = p_phylum, df_phylum = df_phylum)
}

## ---- 7. RESUMEN FINAL -------------------------------------------------------

build_final_summary <- function(perm_int, bd, perm_loc_res) {
  
  cat("\n=== RESUMEN ===\n")
  sink(file.path(OUTDIR, "RESUMEN_estadisticos.txt"))
  cat("=== PERMANOVA Jaccard ~ Species * Locality ===\n"); print(perm_int)
  cat("\n=== Beta dispersion ~ Locality ===\n");           print(bd$bd_test)
  cat("\n=== Beta dispersion ~ Species ===\n");            print(bd$bd_sp_test)
  cat("\n=== PERMANOVA pairwise por localidad ===\n");     print(perm_loc_res)
  sink()
  
  cat("\n====== OBJETIVO 3 COMPLETADO ======\n")
  cat("Outputs en:", file.path(getwd(), OUTDIR), "\n\n")
  cat("Archivos generados:\n")
  cat("  alpha_richness_locality.pdf\n")
  cat("  alpha_diversity_table.csv\n")
  cat("  LM_richness_especie_por_localidad.csv\n")
  cat("  PERMANOVA_Species_x_Locality.txt\n")
  cat("  betadisper_Loc_Sp.txt\n")
  cat("  PERMANOVA_por_localidad_especie.csv\n")
  cat("  PERMANOVA_tabla_unificada.csv\n")
  cat("  PCoA_global_locality_species.pdf\n")
  cat("  PCoA_panel_localidades.pdf\n")
  cat("  Solapamiento_por_rango_taxonomico.pdf / .png / .csv\n")
  cat("  Phylum_by_locality_species.pdf\n")
  cat("  RESUMEN_estadisticos.txt\n")
  
  invisible(NULL)
}


## ---- EXTRA: Porcentajes exactos Ochrophyta y Bacillariophyta ---------------

compute_exact_percentages <- function(df_phylum) {
  
  cat("\n=== Porcentajes Ochrophyta y Bacillariophyta por especie y localidad ===\n")
  
  phyla_interes <- c("Ochrophyta", "Bacillariophyta", "Rhodophyta")
  
  df_exact <- df_phylum |>
    filter(Phylum %in% phyla_interes) |>
    group_by(Sp, Loc, Phylum) |>
    summarise(MeanRel = mean(RelAbund, na.rm = TRUE), .groups = "drop") |>
    mutate(MeanRel_pct = round(100 * MeanRel, 1)) |>
    pivot_wider(names_from = Phylum, values_from = MeanRel_pct, values_fill = 0) |>
    arrange(Sp, Loc)
  
  print(df_exact)
  write.csv(df_exact, 
            file.path(OUTDIR, "Phylum_pct_Ochrophyta_Bacillariophyta.csv"),
            row.names = FALSE)
  
  df_global <- df_phylum |>
    filter(Phylum %in% phyla_interes) |>
    group_by(Sp, Phylum) |>
    summarise(MeanRel_pct = round(100 * mean(RelAbund, na.rm = TRUE), 1),
              SD = round(100 * sd(RelAbund, na.rm = TRUE), 1),
              .groups = "drop")
  
  cat("\nMedia global por especie:\n")
  print(df_global)
  write.csv(df_global,
            file.path(OUTDIR, "Phylum_pct_global_por_especie.csv"),
            row.names = FALSE)
  
  list(by_locality = df_exact, global = df_global)
}

### ══════════════════════════════════════════════════════════════════════════
### Run!
### ══════════════════════════════════════════════════════════════════════════
alpha_df <- build_alpha_df()
plot_alpha_richness(alpha_df)
alpha_loc_res <- run_alpha_lm_by_locality(alpha_df)

beta <- build_beta_matrices()
perm_int <- run_permanova_global(beta)
bd <- run_betadisper_global(beta)
perm_loc_res <- run_permanova_by_locality(beta)
unified_tab <- build_unified_table(perm_int, bd, perm_loc_res)
plot_pcoa_global(beta)

plot_list_pcoa <- list(
  run_pcoa_one_locality(beta, perm_loc_res, "La Bajita"),
  run_pcoa_one_locality(beta, perm_loc_res, "Cercana La Bajita"),
  run_pcoa_one_locality(beta, perm_loc_res, "Charco Verde"),
  run_pcoa_one_locality(beta, perm_loc_res, "El Faro"),
  run_pcoa_one_locality(beta, perm_loc_res, "La Bombilla")
)
names(plot_list_pcoa) <- loc_levels
build_pcoa_locality_panel(plot_list_pcoa)

all_ranks_res <- run_overlap_analysis()
plot_overlap_panel(all_ranks_res)

phylum_result <- plot_phylum_composition()
exact_pct <- compute_exact_percentages(phylum_result$df_phylum)

build_final_summary(perm_int, bd, perm_loc_res)