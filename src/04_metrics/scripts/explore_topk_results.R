#!/usr/bin/env Rscript
# ==============================================================================
# explore_topk_results.R
# Exploratory analysis of top-K rotation invariance results (whole mode).
#
# Complements compute_threshold_free_metrics.R by localizing WHERE rotation
# drift lives in the rankings:
#   1. Coincidence groups (1-4 angles in top-K) vs total gene expression
#   2. Top-K Jaccard by expression decile
#   3. Cross-angle Kendall tau of full score vectors (whole-mode analog of
#      the simulated tau-vs-alpha)
#
# Whole mode only (simulated mode metrics are self-contained).
# Edit SLICES below to control scope.
# Input:  src/01_simulation/outputs/scDesign3/{slice}/whole/data/counts.csv
#         src/03_benchmark/outputs/{slice}/whole/{method}/scdesign3_angle{angle}_results.rds
# Output: src/04_metrics/outputs/{slice}/whole/module1_exploration_topk/{data,figures}/
#         src/04_metrics/outputs/cross_slice/whole/figures/
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# --- Select which slices to run ---
SLICES <- c("anterior1", "anterior2", "posterior1", "posterior2")

TOP_K <- 2000

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# --- Score column and transform per method (nnSVG canonical score: LR_stat,
# aligned with compute_threshold_free_metrics.R) ---
score_config <- list(
  moransi   = list(col = "I",            transform = "identity"),
  spatialde = list(col = "FSV",          transform = "identity"),
  nnsvg     = list(col = "LR_stat",      transform = "identity"),
  sparkx    = list(col = "combinedPval", transform = "neglog10"),
  smash     = list(col = "pval",         transform = "neglog10")
)

# --- Method colors (Wong palette, colorblind-safe) ---
method_colors <- c(
  SPARKX    = "#0072B2",
  NNSVG     = "#D55E00",
  SPATIALDE = "#009E73",
  MORANSI   = "#CC79A7",
  SMASH     = "#F0E442"
)

pair_colors <- c("0\u00b0 vs 30\u00b0" = "#0072B2",
                 "0\u00b0 vs 45\u00b0" = "#D55E00",
                 "0\u00b0 vs 60\u00b0" = "#009E73")

angles <- c(0, 30, 45, 60)
pairs_vs_0 <- list(c("0", "30"), c("0", "45"), c("0", "60"))

# --- Helper: save both PNG (cairo) and PDF ---
save_figure <- function(plot, filename_base, width, height, dpi = 300) {
  ggsave(paste0(filename_base, ".png"), plot, width = width, height = height,
         dpi = dpi, bg = "white")
  ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height, bg = "white")
}

# --- Helper: Jaccard index for two sets ---
jaccard_index <- function(set_a, set_b) {
  inter <- length(intersect(set_a, set_b))
  uni   <- length(union(set_a, set_b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

# --- Containers for cross-slice aggregation ---
all_rotation_tau <- list()

# ==============================================================================
# MAIN LOOP OVER SLICES
# ==============================================================================

for (slice in SLICES) {

  mode <- "whole"

  counts_file <- file.path(project_root, "src", "01_simulation", "outputs",
                           "scDesign3", slice, mode, "data", "counts.csv")
  benchmark_root <- file.path(project_root, "src", "03_benchmark", "outputs",
                              slice, mode)
  expl_dir       <- file.path(project_root, "src", "04_metrics", "outputs",
                              slice, mode, "module1_exploration_topk")
  expl_data_dir  <- file.path(expl_dir, "data")
  expl_fig_dir   <- file.path(expl_dir, "figures")

  if (!file.exists(counts_file)) {
    message("Skipping ", slice, " -- counts.csv not found")
    next
  }
  if (!dir.exists(benchmark_root)) {
    message("Skipping ", slice, " -- benchmark outputs not found")
    next
  }
  for (d in c(expl_data_dir, expl_fig_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  message("\n========== Top-K exploration: ", slice, " / ", mode, " ==========\n")

  # --------------------------------------------------------------------------
  # DATA LOADING
  # --------------------------------------------------------------------------

  # Total counts per gene (harmonized universe: expressed genes only)
  counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
  total_counts <- rowSums(as.matrix(counts))
  gene_df <- data.frame(gene = names(total_counts),
                        total_counts = as.numeric(total_counts),
                        stringsAsFactors = FALSE)
  gene_df <- gene_df[!duplicated(gene_df$gene), ]
  gene_df <- gene_df[gene_df$total_counts > 0, ]

  tools <- list.files(benchmark_root)
  tools <- tools[tools %in% names(score_config)]

  # Scores per tool x angle (full vectors; ranks and bands derived within group)
  score_data <- do.call(rbind, lapply(tools, function(tool) {
    cfg <- score_config[[tool]]
    do.call(rbind, lapply(angles, function(angle) {
      rds_file <- file.path(benchmark_root, tool,
                            paste0("scdesign3_angle", angle, "_results.rds"))
      if (!file.exists(rds_file)) return(NULL)
      df <- readRDS(rds_file)$res_mtest
      if (cfg$transform == "identity") {
        score <- df[[cfg$col]]
      } else if (cfg$transform == "neglog10") {
        score <- -log10(pmax(df[[cfg$col]], .Machine$double.xmin))
      } else {
        stop("Unknown transform '", cfg$transform, "' for tool '", tool, "'")
      }
      data.frame(tool = tool, tool_label = toupper(tool),
                 angle = angle, feature = rownames(df), score = score,
                 stringsAsFactors = FALSE)
    }))
  }))

  if (is.null(score_data) || nrow(score_data) == 0) {
    message("No benchmark data found for ", slice, " -- skipping")
    next
  }

  # Ranks (distinct, via ties.method="first") within tool x angle
  score_data <- score_data %>%
    group_by(tool, angle) %>%
    mutate(rank = rank(-score, ties.method = "first")) %>%
    ungroup()

  # --- Accessors ---
  topk_features <- function(tool, angle) {
    score_data$feature[score_data$tool == tool & score_data$angle == angle &
                         score_data$rank <= TOP_K]
  }
  score_vec <- function(tool, angle) {
    sub <- score_data[score_data$tool == tool & score_data$angle == angle, ]
    v <- sub$score
    names(v) <- sub$feature
    v
  }

  cat("Loaded", nrow(score_data), "score rows for", length(tools),
      "methods x", length(angles), "angles.\n\n")

  # ========================================================================
  # ANALYSIS 1: Coincidence groups (1-4 angles in top-K) vs total counts
  # ========================================================================

  cat("===== ANALYSIS 1: Coincidence vs total counts =====\n\n")

  coincidence_df <- score_data %>%
    filter(rank <= TOP_K) %>%
    distinct(tool, tool_label, feature, angle) %>%
    group_by(tool, tool_label, feature) %>%
    summarise(coincidence = n(), .groups = "drop") %>%
    rename(gene = feature) %>%
    left_join(gene_df, by = "gene")

  write.csv(coincidence_df, file.path(expl_data_dir, "coincidence_counts.csv"),
            row.names = FALSE)

  cat("Top-K union size and coincidence distribution per method:\n")
  print(coincidence_df %>%
          group_by(tool_label, coincidence) %>%
          summarise(n_genes = n(), .groups = "drop") %>%
          tidyr::pivot_wider(names_from = coincidence, values_from = n_genes,
                             names_prefix = "c", values_fill = 0))
  cat("\n")

  coin_plot <- coincidence_df %>%
    mutate(coincidence = factor(coincidence, levels = c("1", "2", "3", "4")))
  p_coincidence <- ggplot(coin_plot,
                          aes(x = coincidence, y = log10(total_counts + 1),
                              fill = coincidence)) +
    geom_violin(alpha = 0.6, scale = "width", trim = TRUE) +
    geom_jitter(width = 0.15, alpha = 0.2, size = 0.8) +
    facet_wrap(~ tool_label, ncol = 3, scales = "free_y") +
    scale_fill_brewer(palette = "Spectral", name = "Coincidences") +
    labs(x = "Number of rotations including the gene in the top-K",
         y = expression(log[10] ~ "total counts"),
         title = paste0("Gene Expression by Top-", TOP_K,
                        " Rotation Coincidence | ", slice),
         subtitle = "Coincidence 4 = stable prioritization | 1-3 = rotation flicker") +
    theme_bw() +
    theme(legend.position = "right",
          strip.text = element_text(face = "bold", size = 11))
  save_figure(p_coincidence, file.path(expl_fig_dir, "coincidence_vs_counts"),
              width = 12, height = 6)

  # ========================================================================
  # ANALYSIS 2: Top-K Jaccard by expression decile
  # ========================================================================

  cat("===== ANALYSIS 2: Jaccard by count decile =====\n\n")

  decile_breaks <- quantile(gene_df$total_counts, probs = seq(0, 1, 0.1))
  decile_breaks <- unique(decile_breaks)
  if (length(decile_breaks) < 11) {
    gene_df$decile <- ntile(gene_df$total_counts, 10)
  } else {
    gene_df$decile <- cut(gene_df$total_counts, breaks = decile_breaks,
                          include.lowest = TRUE, labels = FALSE)
  }
  decile_labels <- gene_df %>%
    group_by(decile) %>%
    summarise(min_c = min(total_counts), max_c = max(total_counts),
              .groups = "drop") %>%
    mutate(label = paste0("D", decile, ": [", round(min_c, 0), ", ",
                          round(max_c, 0), "]"))

  jaccard_decile <- do.call(rbind, lapply(tools, function(tool) {
    do.call(rbind, lapply(pairs_vs_0, function(pair) {
      a <- as.numeric(pair[1]); b <- as.numeric(pair[2])
      do.call(rbind, lapply(sort(unique(gene_df$decile)), function(d) {
        genes_in_decile <- gene_df$gene[gene_df$decile == d]
        sig_a <- intersect(topk_features(tool, a), genes_in_decile)
        sig_b <- intersect(topk_features(tool, b), genes_in_decile)
        data.frame(
          tool = tool, tool_label = toupper(tool), decile = d,
          pair_label = paste0(pair[1], "\u00b0 vs ", pair[2], "\u00b0"),
          jaccard = jaccard_index(sig_a, sig_b),
          n_genes = length(genes_in_decile),
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))

  write.csv(jaccard_decile, file.path(expl_data_dir, "jaccard_by_decile.csv"),
            row.names = FALSE)

  cat("Jaccard by decile (min per method across pairs and deciles D1-D5):\n")
  print(jaccard_decile %>%
          filter(decile <= 5) %>%
          group_by(tool_label) %>%
          summarise(min_j = min(jaccard, na.rm = TRUE), .groups = "drop"))
  cat("\n")

  jaccard_decile <- left_join(jaccard_decile,
                              decile_labels[, c("decile", "label")],
                              by = "decile")
  jaccard_decile$label <- factor(jaccard_decile$label,
                                 levels = decile_labels$label)

  p_jaccard_decile <- ggplot(jaccard_decile,
                             aes(x = label, y = jaccard,
                                 color = pair_label, group = pair_label)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ tool_label, ncol = 3) +
    scale_color_manual(values = pair_colors, name = "Angle pair") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Count decile (total expression threshold)",
         y = paste0("Top-", TOP_K, " Jaccard index"),
         title = paste0("Top-", TOP_K, " Jaccard vs Expression Level | ", slice),
         subtitle = "Rotation consistency of prioritized genes within each expression decile") +
    theme_bw() +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text = element_text(face = "bold", size = 11))
  save_figure(p_jaccard_decile, file.path(expl_fig_dir, "jaccard_by_decile"),
              width = 14, height = 5)

  # ========================================================================
  # ANALYSIS 3: Cross-angle Kendall tau of full score vectors
  # ========================================================================

  cat("===== ANALYSIS 3: Cross-angle Kendall tau =====\n\n")

  rotation_tau <- do.call(rbind, lapply(tools, function(tool) {
    do.call(rbind, lapply(pairs_vs_0, function(pair) {
      a <- as.numeric(pair[1]); b <- as.numeric(pair[2])
      va <- score_vec(tool, a)
      vb <- score_vec(tool, b)[names(va)]
      tau <- suppressWarnings(cor(va, vb, method = "kendall",
                                  use = "complete.obs"))
      data.frame(
        slice = slice, tool = tool, tool_label = toupper(tool),
        angle_a = a, angle_b = b,
        tau_rotation = tau,
        pair_label = paste0(pair[1], "\u00b0 vs ", pair[2], "\u00b0"),
        stringsAsFactors = FALSE
      )
    }))
  }))

  write.csv(rotation_tau, file.path(expl_data_dir, "rotation_tau.csv"),
            row.names = FALSE)

  cat("Cross-angle Kendall tau (0 deg baseline):\n")
  print(tidyr::pivot_wider(rotation_tau,
                           id_cols = tool_label,
                           names_from = pair_label,
                           values_from = tau_rotation))
  cat("\n")

  all_rotation_tau[[length(all_rotation_tau) + 1]] <- rotation_tau

  cat("\nAll exploration outputs saved for ", slice, ".\n")
}

# ==============================================================================
# CROSS-SLICE FLAGSHIP FIGURES
# ==============================================================================

cat("\n===== CROSS-SLICE EXPLORATION FIGURES =====\n\n")

cross_whole_dir <- file.path(project_root, "src", "04_metrics", "outputs",
                             "cross_slice", "whole", "figures")
dir.create(cross_whole_dir, recursive = TRUE, showWarnings = FALSE)

# --- Cross-angle Kendall tau: methods x angle pairs, faceted by slice ---
if (length(all_rotation_tau) > 0) {
  tau_all <- do.call(rbind, all_rotation_tau)
  tau_all$slice <- factor(tau_all$slice, levels = SLICES)

  p_tau_heat <- ggplot(tau_all,
                       aes(x = pair_label, y = tool_label, fill = tau_rotation)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", tau_rotation)),
              size = 3, fontface = "bold") +
    scale_fill_gradient(low = "#F7FCF5", high = "#00441B",
                        name = "Kendall \u03c4") +
    facet_wrap(~ slice, nrow = 1) +
    labs(x = "Angle pair", y = "Method",
         title = "Cross-Angle Kendall \u03c4 of Full Score Vectors (Whole)",
         subtitle = "Ordinal stability of the complete ranking under rotation | Whole-mode analog of the simulated \u03c4 vs \u03b1") +
    theme_bw() +
    theme(strip.text = element_text(face = "bold", size = 11),
          axis.text.x = element_text(angle = 30, hjust = 1))
  save_figure(p_tau_heat, file.path(cross_whole_dir, "rotation_tau_faceted_heatmap"),
              width = 16, height = 5)
  cat("Rotation tau heatmap saved to", cross_whole_dir, "\n")
}

cat("\nAll slices processed.\n")
