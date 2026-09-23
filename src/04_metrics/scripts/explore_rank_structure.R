#!/usr/bin/env Rscript
# ==============================================================================
# explore_rank_structure.R
# Rank-structure exploration answering advisor feedback on the benchmark.
#
#   1. Rank-shift analysis (simulated mode): per-gene change in rank between
#      0 deg and each rotated angle as a function of the alpha signal level.
#      Localizes WHERE in the signal gradient the rotation damage concentrates
#      (e.g. SPARK-X pushes strongest-signal genes down at oblique angles).
#
#   2. Cross-method rank agreement (both modes): 5x5 pairwise Kendall tau
#      between method score vectors at 0 deg. Answers: do different methods
#      produce the same gene ordering under identical conditions?
#
# Note on the Kendall tau used in the main metrics: gene-level tau-b between
# each gene's ranking score and its alpha tag (feature names carry
# "gene_alpha"). Pairs sharing the same alpha are excluded -- never broken by
# an arbitrary rule -- because within-block order carries no ground truth;
# only cross-alpha ordering is scored.
#
# Edit SLICES below to control scope.
# Input:  src/03_benchmark/outputs/{slice}/{mode}/{method}/scdesign3_angle{angle}_results.rds
# Output: src/04_metrics/outputs/{slice}/simulated/module2_rank_structure/{data,figures}/
#         src/04_metrics/outputs/{slice}/{mode}/cross_method_tau/{data,figures}/
# ==============================================================================

library(ggplot2)
library(dplyr)

# --- Select which slices to run ---
SLICES <- c("anterior1", "anterior2", "posterior1", "posterior2")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# --- Score column and transform per method (nnSVG canonical score: LR_stat) ---
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

angles <- c(0, 30, 45, 60)

# --- Helper: save both PNG (cairo) and PDF ---
save_figure <- function(plot, filename_base, width, height, dpi = 300) {
  ggsave(paste0(filename_base, ".png"), plot, width = width, height = height,
         dpi = dpi, bg = "white")
  ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height, bg = "white")
}

# --- Helper: load ranking scores as a named vector (+ alpha tags if simulated) ---
load_scores <- function(benchmark_root, tool, angle, mode) {
  cfg <- score_config[[tool]]
  rds_file <- file.path(benchmark_root, tool,
                        paste0("scdesign3_angle", angle, "_results.rds"))
  if (!file.exists(rds_file)) return(NULL)
  df <- readRDS(rds_file)$res_mtest
  score <- df[[cfg$col]]
  if (cfg$transform == "neglog10") {
    score <- -log10(pmax(score, .Machine$double.xmin))
  }
  out <- data.frame(feature = rownames(df), score = score,
                    stringsAsFactors = FALSE)
  if (mode == "simulated") {
    # Feature names: "gene_alpha"; null-background features ("NULL_gene")
    # parse to alpha = NA and are automatically excluded from the alpha-axis
    # analyses below (band summaries, plots). Their global ranks still
    # participate in the per-tool ranking, which is the correct behavior.
    parts <- strsplit(out$feature, "_", fixed = TRUE)
    out$alpha <- suppressWarnings(as.numeric(sapply(parts, function(x) {
      if (length(x) < 2) return(NA_real_)
      as.numeric(paste(x[-1], collapse = "_"))
    })))
  }
  out
}

# ==============================================================================
# ANALYSIS 1: Rank-shift vs alpha (simulated mode, per slice)
# ==============================================================================

cat("===== ANALYSIS 1: RANK-SHIFT vs ALPHA (simulated) =====\n")

for (slice in SLICES) {

  benchmark_root <- file.path(project_root, "src", "03_benchmark", "outputs",
                              slice, "simulated")
  if (!dir.exists(benchmark_root)) {
    message("Skipping ", slice, "/simulated -- benchmark outputs not found")
    next
  }
  out_dir <- file.path(project_root, "src", "04_metrics", "outputs", slice,
                       "simulated", "module2_rank_structure")
  data_dir <- file.path(out_dir, "data")
  fig_dir  <- file.path(out_dir, "figures")
  for (d in c(data_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  tools <- list.files(benchmark_root)
  tools <- tools[tools %in% names(score_config)]
  if (length(tools) == 0) next

  message("\n--- Rank shift: ", slice, " / simulated ---")

  # Per-gene shift for each rotated angle vs 0 deg, per method
  rank_shift <- do.call(rbind, lapply(tools, function(tool) {
    base <- load_scores(benchmark_root, tool, 0, "simulated")
    if (is.null(base)) return(NULL)
    base$rank_0 <- rank(-base$score, ties.method = "first")
    do.call(rbind, lapply(setdiff(angles, 0), function(angle) {
      rot <- load_scores(benchmark_root, tool, angle, "simulated")
      if (is.null(rot)) return(NULL)
      rot$rank_a <- rank(-rot$score, ties.method = "first")
      m <- merge(base[, c("feature", "alpha", "rank_0")],
                 rot[, c("feature", "rank_a")], by = "feature")
      data.frame(
        slice = slice, tool = tool, tool_label = toupper(tool),
        angle = angle, feature = m$feature, alpha = m$alpha,
        rank_0 = m$rank_0, rank_angle = m$rank_a,
        shift = m$rank_a - m$rank_0,
        stringsAsFactors = FALSE
      )
    }))
  }))

  if (is.null(rank_shift) || nrow(rank_shift) == 0) next
  write.csv(rank_shift, file.path(data_dir, "rank_shift_vs_alpha.csv"),
            row.names = FALSE)

  # Console summary: mean shift per alpha band (0.8-1, 0.5-0.8, 0-0.5) at 45 deg
  cat("\nMean per-gene rank shift at 45 deg by alpha band:\n")
  s45 <- rank_shift %>% filter(angle == 45)
  band_summary <- s45 %>%
    mutate(band = cut(alpha, breaks = c(-0.01, 0.5, 0.8, 1.01),
                      labels = c("[0,0.5]", "(0.5,0.8]", "(0.8,1]"))) %>%
    group_by(tool_label, band) %>%
    summarise(mean_shift = mean(shift), .groups = "drop")
  print(tidyr::pivot_wider(band_summary, id_cols = tool_label,
                           names_from = band, values_from = mean_shift))
  cat("\n")

  # Figure: rank shift vs alpha at 45 deg (worst case), faceted by method.
  # Points = individual genes (gene-level granularity); line = mean shift per
  # alpha level. Positive shift = gene pushed down the ranking at 45 deg.
  d45 <- rank_shift %>% filter(angle == 45)

  p_shift <- ggplot(d45, aes(x = alpha, y = shift)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(color = "grey40", alpha = 0.2, size = 0.7) +
    stat_summary(aes(color = tool_label, group = tool_label),
                 fun = mean, geom = "line", linewidth = 1.3) +
    facet_wrap(~ tool_label, ncol = 5) +
    scale_color_manual(values = method_colors, guide = "none") +
    labs(x = "Signal strength (alpha)", y = "Rank shift (45\u00b0 \u2212 0\u00b0)",
         title = paste0("Where in the Gradient Does Rotation Bite? | ", slice),
         subtitle = "Each point = one gene (alpha-tagged); positive = pushed down at 45\u00b0. Line = mean shift per alpha level. Kendall tau-b scores only cross-alpha ordering (same-alpha pairs are never compared).") +
    theme_bw() +
    theme(strip.text = element_text(face = "bold", size = 11))
  save_figure(p_shift, file.path(fig_dir, "rank_shift_vs_alpha"),
              width = 16, height = 5)

  cat("Rank-shift figure saved for ", slice, ".\n")
}

# ==============================================================================
# ANALYSIS 2: Cross-method rank agreement (both modes, per slice, 0 deg)
# ==============================================================================

cat("\n===== ANALYSIS 2: CROSS-METHOD RANK AGREEMENT (0 deg) =====\n")

for (slice in SLICES) {
  for (mode in c("simulated", "whole")) {

    benchmark_root <- file.path(project_root, "src", "03_benchmark", "outputs",
                                slice, mode)
    if (!dir.exists(benchmark_root)) next

    out_dir <- file.path(project_root, "src", "04_metrics", "outputs", slice,
                         mode, "cross_method_tau")
    data_dir <- file.path(out_dir, "data")
    fig_dir  <- file.path(out_dir, "figures")
    for (d in c(data_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

    tools <- list.files(benchmark_root)
    tools <- tools[tools %in% names(score_config)]
    if (length(tools) < 2) next

    message("\n--- Cross-method tau: ", slice, " / ", mode, " ---")

    # Score vectors at 0 deg (named by feature)
    vecs <- lapply(tools, function(tool) {
      d <- load_scores(benchmark_root, tool, 0, mode)
      if (is.null(d)) return(NULL)
      v <- d$score
      names(v) <- d$feature
      v
    })
    names(vecs) <- tools
    ok <- !sapply(vecs, is.null)
    tools <- tools[ok]
    vecs <- vecs[ok]

    # Pairwise Kendall tau on the shared feature subset
    tau_long <- do.call(rbind, lapply(tools, function(ta) {
      do.call(rbind, lapply(tools, function(tb) {
        common <- intersect(names(vecs[[ta]]), names(vecs[[tb]]))
        tau <- if (ta == tb) 1 else
          suppressWarnings(cor(vecs[[ta]][common], vecs[[tb]][common],
                               method = "kendall"))
        data.frame(slice = slice, mode = mode,
                   tool_a = toupper(ta), tool_b = toupper(tb),
                   tau = tau, n_common = length(common),
                   stringsAsFactors = FALSE)
      }))
    }))

    write.csv(tau_long, file.path(data_dir, "cross_method_tau.csv"),
              row.names = FALSE)

    # Matrix view for printing
    mat <- matrix(tau_long$tau, nrow = length(tools), ncol = length(tools),
                  dimnames = list(toupper(tools), toupper(tools)))
    cat("\nPairwise Kendall tau at 0 deg (", slice, "/", mode, "):\n", sep = "")
    print(round(mat, 3))

    # Heatmap (methods x methods, adaptive green gradient)
    p_cross <- ggplot(tau_long, aes(x = tool_b, y = tool_a, fill = tau)) +
      geom_tile(color = "white", linewidth = 0.8) +
      geom_text(aes(label = sprintf("%.2f", tau)), size = 4, fontface = "bold") +
      scale_fill_gradient(low = "#F7FCF5", high = "#00441B", name = "Kendall tau") +
      labs(x = "Method", y = "Method",
           title = paste0("Cross-Method Rank Agreement at 0\u00b0 | ", slice, " / ", mode),
           subtitle = "Pairwise Kendall tau between full score vectors (shared gene universe per pair)") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    save_figure(p_cross, file.path(fig_dir, "cross_method_tau_heatmap"),
                width = 8, height = 6)

    cat("Cross-method heatmap saved for ", slice, " / ", mode, ".\n")
  }
}

cat("\nAll rank-structure analyses complete.\n")
