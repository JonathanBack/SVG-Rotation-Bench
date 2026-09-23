#!/usr/bin/env Rscript
# ==============================================================================
# explore_architectures.R
# Architecture-aware evaluation of SVG detection (Phase B).
#
# Strategy-2 design principle: cluster for EVALUATION, not input. The raw
# individual genes are fed to the algorithms (unchanged benchmark); the
# architecture clusters only group the results afterwards. Metagenes are used
# purely as a visualization tool, never as algorithm input.
#
# Per slice (simulated mode):
#   1. Cluster the 50 SVGs into spatial architectures: gene-gene Pearson
#      correlation of log1p expression (alpha = 1 features) across spots,
#      Ward.D hierarchical clustering on 1 - correlation distance.
#   2. ComplexHeatmap of the correlation matrix, rows split by cluster.
#   3. Metagene spatial maps: cluster-mean log1p expression on the tissue
#      (visualization only).
#   4. Detection-by-architecture: per method (0 deg), mean rank of each
#      architecture's genes across the weak-signal regime (0 < alpha <= 0.5)
#      -> methods x architectures heatmap (lower mean rank = better
#      detection) + full per-alpha profile CSV and line plot.
#
# Answers: do methods detect certain spatial architectures (e.g. gradients)
# better than others (e.g. hotspots)?
#
# Edit SLICES / K_ARCH below to control scope.
# Input:  src/01_simulation/outputs/scDesign3/{slice}/simulated/data/{counts,location}.csv
#         src/03_benchmark/outputs/{slice}/simulated/{method}/scdesign3_angle0_results.rds
# Output: src/04_metrics/outputs/{slice}/simulated/architecture/{data,figures}/
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# --- Select which slices to run ---
SLICES <- c("anterior1", "anterior2", "posterior1", "posterior2")

# --- Number of spatial-architecture clusters ---
K_ARCH <- 4

# --- Weak-signal regime for the detection-by-architecture aggregation ---
WEAK_ALPHA_MAX <- 0.5

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# --- Score column and transform per method (aligned with the metrics scripts) ---
score_config <- list(
  moransi   = list(col = "I",            transform = "identity"),
  spatialde = list(col = "FSV",          transform = "identity"),
  nnsvg     = list(col = "LR_stat",      transform = "identity"),
  sparkx    = list(col = "combinedPval", transform = "neglog10"),
  smash     = list(col = "pval",         transform = "neglog10")
)

method_colors <- c(
  SPARKX    = "#0072B2",
  NNSVG     = "#D55E00",
  SPATIALDE = "#009E73",
  MORANSI   = "#CC79A7",
  SMASH     = "#F0E442"
)

cluster_colors <- c(
  "A1" = "#0072B2", "A2" = "#D55E00", "A3" = "#009E73", "A4" = "#CC79A7",
  "A5" = "#F0E442", "A6" = "#66C2A5"
)

# --- Helper: save both PNG and PDF ---
save_figure <- function(plot, filename_base, width, height, dpi = 300) {
  ggsave(paste0(filename_base, ".png"), plot, width = width, height = height,
         dpi = dpi, bg = "white")
  ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height, bg = "white")
}

save_heatmap <- function(ht, filename_base, width, height) {
  png(paste0(filename_base, ".png"), width = width, height = height,
      units = "in", res = 300, bg = "white")
  draw(ht, heatmap_legend_side = "bottom")
  dev.off()
  pdf(paste0(filename_base, ".pdf"), width = width, height = height)
  draw(ht, heatmap_legend_side = "bottom")
  dev.off()
}

# ==============================================================================
# MAIN LOOP OVER SLICES
# ==============================================================================

for (slice in SLICES) {

  counts_file <- file.path(project_root, "src", "01_simulation", "outputs",
                           "scDesign3", slice, "simulated", "data", "counts.csv")
  loc_file    <- file.path(project_root, "src", "01_simulation", "outputs",
                           "scDesign3", slice, "simulated", "data", "location.csv")
  benchmark_root <- file.path(project_root, "src", "03_benchmark", "outputs",
                              slice, "simulated")

  if (!file.exists(counts_file) || !dir.exists(benchmark_root)) {
    message("Skipping ", slice, " -- inputs not found")
    next
  }

  out_dir    <- file.path(project_root, "src", "04_metrics", "outputs", slice,
                          "simulated", "architecture")
  data_dir   <- file.path(out_dir, "data")
  fig_dir    <- file.path(out_dir, "figures")
  for (d in c(data_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  message("\n========== Architecture exploration: ", slice, " ==========\n")

  # --------------------------------------------------------------------------
  # STEP 1: Alpha=1 features -> gene-gene correlation -> Ward.D clusters
  # --------------------------------------------------------------------------

  counts <- as.matrix(read.csv(counts_file, row.names = 1, check.names = FALSE))
  fn <- rownames(counts)
  a1_features <- fn[grepl("_1$", fn) & !grepl("^NULL_", fn)]
  if (length(a1_features) < 10) {
    message("Skipping ", slice, " -- too few alpha=1 features")
    next
  }

  expr_a1 <- log1p(counts[a1_features, , drop = FALSE])   # genes x spots
  genes_a1 <- sub("_.*$", "", a1_features)
  rownames(expr_a1) <- genes_a1

  corr_mat <- cor(t(expr_a1))                              # genes x genes

  hc <- hclust(as.dist(1 - corr_mat), method = "ward.D")
  clust <- cutree(hc, k = K_ARCH)
  arch <- paste0("A", clust[genes_a1])
  names(arch) <- genes_a1

  cat("Architecture cluster sizes:\n")
  print(table(arch))
  cat("\n")

  # --------------------------------------------------------------------------
  # STEP 2: Correlation heatmap (ComplexHeatmap, rows split by cluster)
  # --------------------------------------------------------------------------

  arch_factor <- factor(arch[rownames(corr_mat)],
                        levels = paste0("A", 1:K_ARCH))
  col_fun <- colorRamp2(seq(-0.5, 1, length.out = 9),
                        RColorBrewer::brewer.pal(name = "YlGnBu", n = 9))

  # NOTE: row_split as a factor with cluster_rows = TRUE (dendrograms are
  # drawn within each architecture slice); passing the external hclust object
  # together with a factor split is not supported by ComplexHeatmap.
  ht <- Heatmap(
    corr_mat,
    column_title = paste0(slice, ": gene-gene correlation (alpha = 1)"),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_split = arch_factor,
    col = col_fun,
    row_gap = unit(0.5, "mm"),
    column_gap = unit(0.5, "mm"),
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_parent_dend_line = FALSE,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    heatmap_legend_param = list(title = "Correlation",
                                title_position = "topcenter",
                                direction = "horizontal")
  )
  save_heatmap(ht, file.path(fig_dir, "architecture_correlation_heatmap"),
               width = 7, height = 6)

  # --------------------------------------------------------------------------
  # STEP 3: Metagene spatial maps (visualization only, never algorithm input)
  # --------------------------------------------------------------------------

  loc <- read.csv(loc_file, row.names = 1, check.names = FALSE)
  spots <- colnames(counts)
  loc <- loc[spots, ]

  metagene_df <- do.call(rbind, lapply(levels(arch_factor), function(a) {
    g <- names(arch)[arch == a]
    if (length(g) == 0) return(NULL)
    m <- colMeans(log1p(counts[paste0(g, "_1"), , drop = FALSE]))
    data.frame(spot = spots, spatial1 = loc$spatial1, spatial2 = loc$spatial2,
               metagene = m, architecture = a, n_genes = length(g),
               stringsAsFactors = FALSE)
  }))

  write.csv(metagene_df, file.path(data_dir, "metagene_maps.csv"), row.names = FALSE)

  p_meta <- ggplot(metagene_df, aes(x = spatial1, y = spatial2, color = metagene)) +
    geom_point(size = 0.9) +
    scale_colour_gradientn(colors = RColorBrewer::brewer.pal(name = "YlGnBu", n = 9)) +
    facet_wrap(~ architecture) +
    labs(x = NULL, y = NULL,
         title = paste0("Spatial Architectures (metagenes, alpha = 1) | ", slice),
         subtitle = "Cluster-mean log1p expression -- visualization only, never fed to the algorithms") +
    theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          strip.text = element_text(face = "bold", size = 11))
  save_figure(p_meta, file.path(fig_dir, "metagene_maps"), width = 10, height = 8)

  # --------------------------------------------------------------------------
  # STEP 4: Detection-by-architecture (methods x architectures, 0 deg)
  # --------------------------------------------------------------------------

  tools <- list.files(benchmark_root)
  tools <- tools[tools %in% names(score_config)]

  rank_profiles <- do.call(rbind, lapply(tools, function(tool) {
    cfg <- score_config[[tool]]
    rds_file <- file.path(benchmark_root, tool, "scdesign3_angle0_results.rds")
    if (!file.exists(rds_file)) return(NULL)
    df <- readRDS(rds_file)$res_mtest
    score <- if (cfg$transform == "identity") df[[cfg$col]] else
      -log10(pmax(df[[cfg$col]], .Machine$double.xmin))

    feat <- rownames(df)
    is_null <- grepl("^NULL_", feat)
    gene <- sub("_.*$", "", feat)
    parts <- strsplit(feat, "_", fixed = TRUE)
    alpha <- suppressWarnings(as.numeric(sapply(parts, function(x) {
      if (length(x) < 2) return(NA_real_)
      as.numeric(paste(x[-1], collapse = "_"))
    })))
    alpha[is_null] <- NA_real_

    keep <- !is_null & !is.na(alpha) & gene %in% names(arch)
    rk <- rank(-score, ties.method = "first")

    data.frame(
      tool = tool, tool_label = toupper(tool),
      gene = gene[keep], alpha = alpha[keep],
      architecture = arch[gene[keep]],
      rank = rk[keep],
      n_features = length(feat),
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(rank_profiles) || nrow(rank_profiles) == 0) {
    message("No benchmark profiles for ", slice, " -- skipping detection analysis")
    next
  }

  # Per (tool, architecture, alpha): mean rank of the cluster's genes
  profile_agg <- rank_profiles %>%
    group_by(tool_label, architecture, alpha) %>%
    summarise(mean_rank = mean(rank), n = n(), .groups = "drop")
  write.csv(profile_agg, file.path(data_dir, "detection_by_architecture.csv"),
            row.names = FALSE)

  # Weak-signal regime aggregation (0 < alpha <= 0.5)
  weak <- profile_agg %>% filter(alpha > 0, alpha <= WEAK_ALPHA_MAX) %>%
    group_by(tool_label, architecture) %>%
    summarise(mean_rank = mean(mean_rank), .groups = "drop")

  cat("Detection-by-architecture (mean rank, 0 < alpha <= 0.5; lower = better):\n")
  print(tidyr::pivot_wider(weak, id_cols = tool_label,
                           names_from = architecture, values_from = mean_rank))
  cat("\n")

  # Heatmap: methods x architectures. Lower mean rank = better detection, so
  # the green gradient is reversed (dark green = low rank = good).
  p_detect <- ggplot(weak, aes(x = architecture, y = tool_label, fill = mean_rank)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.0f", mean_rank)), size = 4, fontface = "bold") +
    scale_fill_gradient(low = "#00441B", high = "#F7FCF5", name = "Mean rank") +
    labs(x = "Spatial architecture", y = "Method",
         title = paste0("Detection by Spatial Architecture | ", slice),
         subtitle = paste0("Mean rank of each architecture's genes, weak-signal regime (0 < alpha <= ",
                           WEAK_ALPHA_MAX, ") | lower = better detection")) +
    theme_bw() +
    theme(axis.text.x = element_text(face = "bold"))
  save_figure(p_detect, file.path(fig_dir, "detection_by_architecture_heatmap"),
              width = 8, height = 5)

  # Line plot: full per-alpha profile
  p_lines <- ggplot(profile_agg %>% filter(alpha > 0),
                    aes(x = alpha, y = mean_rank, color = architecture,
                        group = architecture)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    facet_wrap(~ tool_label, ncol = 3) +
    scale_color_manual(values = cluster_colors, name = "Architecture") +
    labs(x = "Signal strength (alpha)", y = "Mean rank",
         title = paste0("Architecture Detection Profile Across Signal Levels | ", slice),
         subtitle = "Mean rank of each architecture's genes per alpha level (lower = better)") +
    theme_bw() +
    theme(strip.text = element_text(face = "bold", size = 11))
  save_figure(p_lines, file.path(fig_dir, "detection_by_architecture_lines"),
              width = 12, height = 6)

  cat("Architecture outputs saved for ", slice, ".\n")
}

cat("\nAll architecture analyses complete.\n")
