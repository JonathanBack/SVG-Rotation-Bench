#!/usr/bin/env Rscript
# ==============================================================================
# explore_module1_results.R
# Exploratory analysis of Module 1 (rotation invariance) results for the whole
# gene set. Links rotation-consistency patterns to gene expression level.
#
# Two analyses per slice:
#   1. Coincidence groups (1-4) vs total gene expression
#      For each method, genes are classified by how many angles called them
#      significant (1 = one angle only, 4 = all four angles). Violin plot of
#      log10(total counts) by coincidence group, faceted by method.
#   2. Jaccard index by count decile
#      Genes binned into 10 deciles by total expression. Jaccard computed
#      per-decile to show whether high-expression genes are more consistent
#      across rotations. Decile labels show actual count thresholds.
#
# Whole mode only (simulated mode has alpha as its own signal dimension).
# Per-slice only (no cross-slice aggregation yet).
#
# Edit SLICES below to control scope.
# Input:  src/01_simulation/outputs/scDesign3/{slice}/whole/data/counts.csv
#         src/03_benchmark/outputs/{slice}/whole/{method}/scdesign3_angle{angle}_results.rds
# Output: src/04_metrics/outputs/{slice}/whole/module1_exploration/{data,figures}/
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# --- Select which slices to run ---
SLICES <- c("anterior1", "anterior2", "posterior1", "posterior2")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# --- Helpers ---
save_figure <- function(plot, filename_base, width, height, dpi = 300) {
  ggsave(paste0(filename_base, ".png"), plot, width = width, height = height,
         dpi = dpi, bg = "white")
  ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height, bg = "white")
}

jaccard_index <- function(set_a, set_b) {
  inter <- length(intersect(set_a, set_b))
  uni   <- length(union(set_a, set_b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

# --- Adjusted p-value columns per method ---
adj_pval_col <- list(
  sparkx    = "adjustedPval",
  nnsvg     = "padj",
  spatialde = "qval",
  moransi   = "pval_sim_fdr_bh",
  smash     = "adjusted_pval"
)

angles <- c(0, 30, 45, 60)
pairs_vs_0 <- list(c("0", "30"), c("0", "45"), c("0", "60"))

for (slice in SLICES) {

  mode <- "whole"

  message("\n========== Exploring Module 1: ", slice, " / ", mode, " ==========\n")

  # Paths
  counts_file <- file.path(project_root, "src", "01_simulation", "outputs",
                           "scDesign3", slice, mode, "data", "counts.csv")
  benchmark_root <- file.path(project_root, "src", "03_benchmark", "outputs",
                              slice, mode)
  expl_dir       <- file.path(project_root, "src", "04_metrics", "outputs",
                              slice, mode, "module1_exploration")
  expl_data_dir  <- file.path(expl_dir, "data")
  expl_fig_dir   <- file.path(expl_dir, "figures")
  dir.create(expl_data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(expl_fig_dir,  recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(counts_file)) {
    message("Skipping ", slice, " -- counts.csv not found")
    next
  }
  if (!dir.exists(benchmark_root)) {
    message("Skipping ", slice, " -- benchmark outputs not found")
    next
  }

  # Load counts matrix (genes x spots) and compute total expression per gene
  counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
  total_counts <- rowSums(counts)
  gene_df <- data.frame(
    gene = names(total_counts),
    total_counts = as.numeric(total_counts),
    stringsAsFactors = FALSE
  )

  # Identify available methods
  tools <- list.files(benchmark_root)
  tools <- tools[tools %in% names(adj_pval_col)]

  if (length(tools) == 0) {
    message("No benchmark methods found for ", slice, " -- skipping")
    next
  }

  # ========================================================================
  # ANALYSIS 1: Coincidence groups (1-4) vs total counts
  # ========================================================================

  coincidence_df <- do.call(rbind, lapply(tools, function(tool) {
    pcol <- adj_pval_col[[tool]]

    # Gather significant genes per angle
    sig_by_angle <- lapply(angles, function(angle) {
      rds_file <- file.path(benchmark_root, tool,
                            paste0("scdesign3_angle", angle, "_results.rds"))
      if (!file.exists(rds_file)) return(character(0))
      df <- readRDS(rds_file)$res_mtest
      adj_p <- df[[pcol]]
      rownames(df)[!is.na(adj_p) & adj_p < 0.05]
    })
    names(sig_by_angle) <- as.character(angles)

    # Count how many angles called each gene significant (0-4)
    n_sig <- sapply(gene_df$gene, function(g) {
      sum(sapply(sig_by_angle, function(s) g %in% s))
    })

    data.frame(
      tool = tool,
      gene = gene_df$gene,
      total_counts = gene_df$total_counts,
      coincidence = n_sig,
      stringsAsFactors = FALSE
    )
  }))

  # Exclude 0-coincidence group (never significant in any angle)
  coincidence_df <- coincidence_df %>% filter(coincidence > 0)
  coincidence_df$coincidence <- factor(coincidence_df$coincidence,
                                       levels = c("1", "2", "3", "4"))
  coincidence_df$tool_label <- toupper(coincidence_df$tool)

  write.csv(coincidence_df,
            file.path(expl_data_dir, "coincidence_counts.csv"),
            row.names = FALSE)

  # Violin + jitter plot (log10 with +1 pseudocount to handle zeros)
  p_coincidence <- ggplot(coincidence_df,
                          aes(x = coincidence, y = log10(total_counts + 1),
                              fill = coincidence)) +
    geom_violin(alpha = 0.6, scale = "width", trim = TRUE) +
    geom_jitter(width = 0.15, alpha = 0.2, size = 0.8) +
    facet_wrap(~ tool_label, ncol = 3, scales = "free_y") +
    scale_fill_brewer(palette = "Spectral", name = "Coincidences") +
    labs(x = "Number of angles calling gene significant",
         y = expression(log[10] ~ "total counts"),
         title = paste0("Gene expression by rotation-consistency | ", slice),
         subtitle = "Higher coincidence = detected across more rotations") +
    theme_bw() +
    theme(legend.position = "right",
          strip.text = element_text(face = "bold", size = 11))

  save_figure(p_coincidence,
              file.path(expl_fig_dir, "coincidence_vs_counts"),
              width = 12, height = 6)

  cat("Coincidence plot saved.\n")

  # ========================================================================
  # ANALYSIS 2: Jaccard by count decile
  # ========================================================================

  # Compute decile thresholds once per slice (shared across methods)
  decile_breaks <- quantile(gene_df$total_counts, probs = seq(0, 1, 0.1))
  # Remove duplicates to avoid empty bins
  decile_breaks <- unique(decile_breaks)
  # If fewer than 11 unique breaks, fall back to ntile from dplyr
  if (length(decile_breaks) < 11) {
    gene_df$decile <- ntile(gene_df$total_counts, 10)
  } else {
    gene_df$decile <- cut(gene_df$total_counts,
                          breaks = decile_breaks,
                          include.lowest = TRUE,
                          labels = FALSE)
  }

  # Build decile labels showing actual thresholds
  decile_labels <- gene_df %>%
    group_by(decile) %>%
    summarise(
      min_c = min(total_counts),
      max_c = max(total_counts),
      .groups = "drop"
    ) %>%
    mutate(label = paste0(decile, ": [", round(min_c, 0), ", ", round(max_c, 0), "]"))

  jaccard_decile <- do.call(rbind, lapply(tools, function(tool) {
    pcol <- adj_pval_col[[tool]]

    sig_by_angle <- lapply(angles, function(angle) {
      rds_file <- file.path(benchmark_root, tool,
                            paste0("scdesign3_angle", angle, "_results.rds"))
      if (!file.exists(rds_file)) return(character(0))
      df <- readRDS(rds_file)$res_mtest
      adj_p <- df[[pcol]]
      rownames(df)[!is.na(adj_p) & adj_p < 0.05]
    })
    names(sig_by_angle) <- as.character(angles)

    do.call(rbind, lapply(pairs_vs_0, function(pair) {
      a <- pair[1]; b <- pair[2]

      do.call(rbind, lapply(sort(unique(gene_df$decile)), function(d) {
        genes_in_decile <- gene_df$gene[gene_df$decile == d]
        sig_a <- intersect(sig_by_angle[[a]], genes_in_decile)
        sig_b <- intersect(sig_by_angle[[b]], genes_in_decile)

        data.frame(
          tool = tool, tool_label = toupper(tool),
          decile = d,
          pair_label = paste0(a, "\u00b0 vs ", b, "\u00b0"),
          jaccard = jaccard_index(sig_a, sig_b),
          n_genes = length(genes_in_decile),
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))

  # Attach decile threshold labels
  jaccard_decile <- left_join(
    jaccard_decile,
    decile_labels[, c("decile", "label")],
    by = "decile"
  )
  # Order deciles naturally
  jaccard_decile$label <- factor(jaccard_decile$label,
                                  levels = decile_labels$label)

  write.csv(jaccard_decile,
            file.path(expl_data_dir, "jaccard_by_decile.csv"),
            row.names = FALSE)

  pair_colors <- c("0\u00b0 vs 30\u00b0" = "#0072B2",
                   "0\u00b0 vs 45\u00b0" = "#D55E00",
                   "0\u00b0 vs 60\u00b0" = "#009E73")

  p_jaccard_decile <- ggplot(jaccard_decile,
                             aes(x = label, y = jaccard,
                                 color = pair_label, group = pair_label)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ tool_label, ncol = 3) +
    scale_color_manual(values = pair_colors, name = "Angle pair") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Count decile (total expression threshold)",
         y = "Jaccard index",
         title = paste0("Jaccard vs Expression Level | ", slice),
         subtitle = "Rotation consistency within each expression decile") +
    theme_bw() +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text = element_text(face = "bold", size = 11))

  save_figure(p_jaccard_decile,
              file.path(expl_fig_dir, "jaccard_by_decile"),
              width = 14, height = 5)

  cat("Jaccard-by-decile plot saved.\n")
  cat("\nExploration outputs saved for ", slice, ".\n")
}

cat("\nAll slices processed.\n")