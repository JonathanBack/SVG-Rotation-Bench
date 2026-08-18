#!/usr/bin/env Rscript
# ==============================================================================
# compute_all_metrics.R
# Metrics computation for SVG rotation invariance benchmark.
#
# Kept metrics (per slice):
#   - Jaccard index + Jaccard heatmap (all angles)            [both modes]
#   - Venn diagrams across 4 angles                           [both modes]
#   - FPR at alpha = 0, 0 deg only                            [simulated only]
#   - Sensitivity by alpha, 0 deg only                        [simulated only]
#   - Confusion matrix heatmap at 0 deg                       [simulated only]
#   - Runtime barplot (all angles)                            [both modes]
#
# Per-slice only (no cross-slice aggregation). Module 2 is skipped
# entirely when mode == "whole" because no alpha ground truth is available.
#
# Edit SLICES and MODES below to control scope.
# Input:  src/03_benchmark/outputs/{slice}/{mode}/{method}/scdesign3_angle{angle}_results.rds
# Output: src/04_metrics/outputs/{slice}/{mode}/{module1_rotation|module2_performance}/...
# ==============================================================================

# ==============================================================================
# SECTION 0: SETUP AND CONFIGURATION
# ==============================================================================

library(ggplot2)
library(reshape2)
library(dplyr)
library(ggvenn)

# --- Select which slices and modes to run ---
SLICES <- c("anterior1", "anterior2", "posterior1", "posterior2")
MODES  <- c("simulated", "whole")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# --- Score configuration: which column to use as the ranking "statistic" ---
# negate = TRUE means we flip the sign so that higher score = more significant
score_config <- list(
  sparkx    = list(col = "combinedPval", negate = TRUE),
  nnsvg     = list(col = "LR_stat",      negate = FALSE),
  spatialde = list(col = "FSV",          negate = FALSE),
  moransi   = list(col = "I",            negate = FALSE),
  smash     = list(col = "pval",         negate = TRUE)
)

# --- Adjusted p-value column for significance filtering (adj p < 0.05) ---
adj_pval_col <- list(
  sparkx    = "adjustedPval",
  nnsvg     = "padj",
  spatialde = "qval",
  moransi   = "pval_sim_fdr_bh",
  smash     = "adjusted_pval"
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
pairs_all <- list(c("0", "30"), c("0", "45"), c("0", "60"),
                  c("30", "45"), c("30", "60"), c("45", "60"))
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

# ==============================================================================
# MAIN LOOP OVER SLICES AND MODES
# ==============================================================================

for (slice in SLICES) {
  for (mode in MODES) {

    benchmark_root <- file.path(project_root, "src", "03_benchmark", "outputs",
                                slice, mode)
    metrics_root   <- file.path(project_root, "src", "04_metrics", "outputs",
                                slice, mode)

    if (!dir.exists(benchmark_root)) {
      message("Skipping ", slice, "/", mode, " -- benchmark outputs not found")
      next
    }

    module1_dir       <- file.path(metrics_root, "module1_rotation")
    module1_data_dir  <- file.path(module1_dir, "data")
    module1_fig_dir   <- file.path(module1_dir, "figures")
    module2_dir       <- file.path(metrics_root, "module2_performance")
    module2_data_dir  <- file.path(module2_dir, "data")
    module2_fig_dir   <- file.path(module2_dir, "figures")

    for (d in c(module1_data_dir, module1_fig_dir,
                module2_data_dir, module2_fig_dir)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }

    message("\n========== Metrics: ", slice, " / ", mode, " ==========\n")

    # --------------------------------------------------------------------------
    # SECTION 1: DATA LOADING
    # --------------------------------------------------------------------------

    cat("Loading benchmark results for", slice, "/", mode, "...\n")

    tools <- list.files(benchmark_root)
    tools <- tools[tools %in% names(score_config)]

    all_data <- do.call(rbind, lapply(tools, function(tool) {
      cfg  <- score_config[[tool]]
      pcol <- adj_pval_col[[tool]]
      if (is.null(cfg) || is.null(pcol)) return(NULL)

      do.call(rbind, lapply(angles, function(angle) {
        rds_file <- file.path(benchmark_root, tool,
                              paste0("scdesign3_angle", angle, "_results.rds"))
        runtime_file <- file.path(benchmark_root, tool,
                                  paste0("scdesign3_angle", angle, "_runtime.csv"))
        if (!file.exists(rds_file)) return(NULL)

        df <- readRDS(rds_file)$res_mtest
        feature_names <- rownames(df)

        # Feature names follow "gene_alpha" for simulated mode; plain gene
        # symbols for whole mode (no alpha suffix => NA).
        parts <- strsplit(feature_names, "_")
        gene  <- sapply(parts, `[`, 1)
        alpha <- suppressWarnings(as.numeric(sapply(parts, function(x) {
          if (length(x) < 2) return(NA_real_)
          as.numeric(paste(x[-1], collapse = "_"))
        })))

        raw_score <- df[[cfg$col]]
        score <- if (cfg$negate) -raw_score else raw_score
        adj_p <- df[[pcol]]

        runtime_val <- if (file.exists(runtime_file)) read.csv(runtime_file)$elapsed_sec else NA_real_

        data.frame(
          tool = tool, angle = angle, gene = gene, alpha = alpha,
          feature = feature_names, score = score, adj_p = adj_p,
          significant = !is.na(adj_p) & adj_p < 0.05,
          runtime_sec = runtime_val,
          stringsAsFactors = FALSE
        )
      }))
    }))

    if (is.null(all_data) || nrow(all_data) == 0) {
      message("No benchmark data found for ", slice, "/", mode, " -- skipping")
      next
    }

    all_data$tool_label <- toupper(all_data$tool)
    all_data$angle <- as.numeric(all_data$angle)

    cat("Loaded", nrow(all_data), "rows for",
        length(unique(all_data$tool)), "methods across",
        length(unique(all_data$angle)), "angles.\n")
    cat("Methods:", paste(unique(all_data$tool_label), collapse = ", "), "\n")
    if (mode == "simulated") {
      cat("Alpha levels:", paste(sort(unique(all_data$alpha)), collapse = ", "), "\n")
    }
    cat("\n")

    # ==============================================================================
    # MODULE 1: ROTATION INVARIANCE  (both modes)
    # ==============================================================================
    # Goal: Evaluate if methods yield the same result across rotated coordinates
    # (0, 30, 45, 60 degrees). Uses Jaccard index and Venn diagrams over the
    # significant gene sets.
    # ==============================================================================

    cat("===== MODULE 1: ROTATION INVARIANCE =====\n\n")

    # --- Metric 1.1: Set Overlap (Jaccard & Venn) ---
    sig_sets_by_tool <- lapply(tools, function(tool) {
      sets <- lapply(angles, function(angle) {
        subset <- all_data[all_data$tool == tool & all_data$angle == angle, ]
        subset$feature[subset$significant]
      })
      names(sets) <- as.character(angles)
      sets
    })
    names(sig_sets_by_tool) <- tools

    jaccard_results <- do.call(rbind, lapply(tools, function(tool) {
      sets <- sig_sets_by_tool[[tool]]
      do.call(rbind, lapply(pairs_all, function(pair) {
        a <- pair[1]; b <- pair[2]
        if (!a %in% names(sets) || !b %in% names(sets)) return(NULL)
        data.frame(
          tool = tool, tool_label = toupper(tool),
          angle_a = as.numeric(a), angle_b = as.numeric(b),
          jaccard = jaccard_index(sets[[a]], sets[[b]]),
          pair_label = paste0(a, "\u00b0 vs ", b, "\u00b0"),
          stringsAsFactors = FALSE
        )
      }))
    }))
    if (!is.null(jaccard_results) && nrow(jaccard_results) > 0) {
      write.csv(jaccard_results, file.path(module1_data_dir, "jaccard.csv"), row.names = FALSE)
    }

    # Venn diagrams (4-way: 0, 30, 45, 60) per method
    # Per-set semi-transparent fills (Wong palette), black borders,
    # labels show both count and percent (1 decimal digit).

    for (tool in tools) {
      sets <- sig_sets_by_tool[[tool]]
      if (length(sets) < 4) next
      venn_data <- list(
        Original   = sets[["0"]],
        `30 degree` = sets[["30"]],
        `45 degree` = sets[["45"]],
        `60 degree` = sets[["60"]]
      )
      p_venn <- ggvenn(
        venn_data,
        fill_color = c("#0072B2", "#D55E00", "#009E73", "#CC79A7"),
        fill_alpha = 0.35,
        stroke_color = "black",
        stroke_size  = 0.8,
        show_stats = "cp",
        digits = 1,
        text_color = "black"
      ) +
        labs(title = paste0(toupper(tool), ": Significant SVG Overlap Across Rotations (", slice, "/", mode, ")"),
             subtitle = "Adjusted p < 0.05 | Intersection = consistent across all angles") +
        theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
      save_figure(p_venn, file.path(module1_fig_dir, paste0("venn_", tool)),
                  width = 10, height = 8)
    }

    # Cross-method Jaccard heatmap
    if (!is.null(jaccard_results) && nrow(jaccard_results) > 0) {
      p_jaccard <- ggplot(jaccard_results,
                          aes(x = pair_label, y = tool_label, fill = jaccard)) +
        geom_tile(color = "white", linewidth = 0.8) +
        geom_text(aes(label = sprintf("%.3f", jaccard)), size = 3.5, fontface = "bold") +
        scale_fill_gradient2(low = "#D73027", mid = "#FFFFCC", high = "#1A9850",
                             midpoint = 0.95, limits = c(0, 1), name = "Jaccard index") +
        labs(x = "Angle pair", y = "Method",
             title = paste0("Set Overlap: Jaccard Index Across Angle Pairs (", slice, "/", mode, ")"),
             subtitle = "Jaccard = |intersection| / |union| of significant gene sets") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
      save_figure(p_jaccard, file.path(module1_fig_dir, "jaccard_heatmap"),
                  width = 10, height = 6)
    }

    cat("Module 1 outputs saved to", module1_dir, "\n\n")

    # ==============================================================================
    # MODULE 2: STATISTICAL PERFORMANCE  (simulated mode only)
    # ==============================================================================
    # Ground truth: alpha = 0 is pure noise (TN); alpha > 0 is spatial signal (TP).
    # Only evaluated at the 0 deg baseline, since invariance is already covered
    # by Module 1 and alpha-ground-truth requires the simulation.
    # ==============================================================================

    if (mode != "simulated") {
      cat("Skipping Module 2 (no alpha ground truth for mode = ", mode, ")\n\n")
    } else {
      cat("===== MODULE 2: STATISTICAL PERFORMANCE (0 deg only) =====\n\n")

      baseline_data <- all_data %>% filter(angle == 0)

      # --- Confusion matrix counts at adj p < 0.05 ---
      class_df <- do.call(rbind, lapply(tools, function(tool) {
        sub <- all_data[all_data$tool == tool & all_data$angle == 0, ]
        if (nrow(sub) == 0) return(NULL)
        truth <- as.integer(sub$alpha > 0)
        predicted <- as.integer(sub$significant)
        data.frame(
          tool = tool, tool_label = toupper(tool), angle = 0,
          TP = sum(truth == 1 & predicted == 1, na.rm = TRUE),
          TN = sum(truth == 0 & predicted == 0, na.rm = TRUE),
          FP = sum(truth == 0 & predicted == 1, na.rm = TRUE),
          FN = sum(truth == 1 & predicted == 0, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }))

      # --- Metric 2.1: Statistical Calibration (FPR) ---
      fpr_summary <- baseline_data %>%
        filter(alpha == 0) %>%
        group_by(tool, tool_label) %>%
        summarise(
          n_negatives = n(),
          n_false_pos = sum(significant, na.rm = TRUE),
          fpr = mean(significant, na.rm = TRUE),
          .groups = "drop"
        )
      write.csv(fpr_summary, file.path(module2_data_dir, "fpr.csv"), row.names = FALSE)

      cat("FPR (alpha = 0, angle = 0):\n")
      print(fpr_summary)

      p_fpr <- ggplot(fpr_summary, aes(x = tool_label, y = fpr, fill = tool_label)) +
        geom_col(color = "black", linewidth = 0.3) +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
        geom_text(aes(label = sprintf("%.3f", fpr)), vjust = -0.5, size = 4, fontface = "bold") +
        scale_fill_manual(values = method_colors, guide = "none") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        labs(x = "Method", y = "False Positive Rate",
             title = paste0("Statistical Calibration: FPR at alpha = 0 (0\u00b0 baseline) | ", slice),
             subtitle = "Dashed line = nominal 0.05 threshold | FPR >> 0.05 = overly aggressive") +
        theme_bw()
      save_figure(p_fpr, file.path(module2_fig_dir, "fpr_barplot"), width = 9, height = 6)

      # --- Metric 2.2: Limit of Detection (Sensitivity across signal gradients)
      sensitivity_by_alpha <- baseline_data %>%
        filter(alpha > 0) %>%
        group_by(tool, tool_label, alpha) %>%
        summarise(
          n_total = n(),
          n_significant = sum(significant, na.rm = TRUE),
          sensitivity = mean(significant, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(tool, alpha)
      write.csv(sensitivity_by_alpha,
                file.path(module2_data_dir, "sensitivity_by_alpha.csv"), row.names = FALSE)

      cat("\nSensitivity by alpha (angle = 0):\n")
      print(sensitivity_by_alpha)

      p_sens <- ggplot(sensitivity_by_alpha,
                       aes(x = alpha, y = sensitivity, color = tool_label, group = tool_label)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 2.5) +
        scale_color_manual(values = method_colors, name = "Method") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                           limits = c(0, 1.05)) +
        labs(x = "Signal strength (alpha)", y = "Sensitivity (Power)",
             title = paste0("Limit of Detection: Sensitivity Across Signal Gradients (0\u00b0) | ", slice),
             subtitle = "At what signal strength does each method recognize an SVG?") +
        theme_bw() +
        theme(legend.position = "right")
      save_figure(p_sens, file.path(module2_fig_dir, "sensitivity_lines"), width = 10, height = 6)

      # --- Confusion Matrix at 0 deg ---
      if (!is.null(class_df) && nrow(class_df) > 0) {
        conf_melt <- melt(class_df, id.vars = c("tool_label", "angle"),
                          measure.vars = c("TP", "TN", "FP", "FN"),
                          variable.name = "Cell", value.name = "Count")

        p_conf <- ggplot(conf_melt, aes(x = Cell, y = reorder(tool_label, Count), fill = Count)) +
          geom_tile(color = "white", linewidth = 1) +
          geom_text(aes(label = Count), size = 5, fontface = "bold") +
          scale_fill_gradient(low = "#FEE8C8", high = "#E6550D", name = "Count") +
          labs(x = "Confusion Matrix Cell", y = "Method",
               title = paste0("Confusion Matrix at 0\u00b0 - All Methods | ", slice),
               subtitle = "Adjusted p < 0.05 significance threshold") +
          theme_bw() +
          theme(axis.text.y = element_text(face = "bold", size = 11))
        save_figure(p_conf, file.path(module2_fig_dir, "confusion_matrix_heatmap"),
                    width = 10, height = 6)
      }

      # --- Metric 2.4: Jaccard by alpha (rotation consistency vs signal strength) ---
      # For each alpha level, compute Jaccard between 0° and each rotated angle
      # using ONLY features at that alpha. Shows how signal strength affects
      # cross-angle set agreement.

      jaccard_by_alpha <- do.call(rbind, lapply(tools, function(tool) {
        do.call(rbind, lapply(pairs_vs_0, function(pair) {
          a <- pair[1]; b <- pair[2]
          do.call(rbind, lapply(sort(unique(all_data$alpha)), function(alpha_level) {
            sig_a <- all_data$feature[
              all_data$tool == tool & all_data$angle == as.numeric(a) &
              all_data$alpha == alpha_level & all_data$significant
            ]
            sig_b <- all_data$feature[
              all_data$tool == tool & all_data$angle == as.numeric(b) &
              all_data$alpha == alpha_level & all_data$significant
            ]
            data.frame(
              tool = tool, tool_label = toupper(tool),
              alpha = alpha_level, angle_a = a, angle_b = b,
              jaccard = jaccard_index(sig_a, sig_b),
              pair_label = paste0(a, "\u00b0 vs ", b, "\u00b0"),
              stringsAsFactors = FALSE
            )
          }))
        }))
      }))
      write.csv(jaccard_by_alpha,
                file.path(module2_data_dir, "jaccard_by_alpha.csv"), row.names = FALSE)

      cat("\nJaccard by alpha (0\u00b0 vs rotated):\n")
      print(jaccard_by_alpha %>%
              filter(!is.na(jaccard)) %>%
              group_by(tool_label) %>%
              summarise(min_j = min(jaccard), max_j = max(jaccard), mean_j = mean(jaccard),
                        .groups = "drop"))

      pair_colors <- c("0\u00b0 vs 30\u00b0" = "#0072B2",
                       "0\u00b0 vs 45\u00b0" = "#D55E00",
                       "0\u00b0 vs 60\u00b0" = "#009E73")

      p_jaccard_alpha <- ggplot(jaccard_by_alpha,
                                aes(x = alpha, y = jaccard,
                                    color = pair_label, group = pair_label)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        facet_wrap(~ tool_label, ncol = 3) +
        scale_color_manual(values = pair_colors, name = "Angle pair") +
        scale_y_continuous(limits = c(0, 1)) +
        labs(x = "Signal strength (alpha)", y = "Jaccard index",
             title = paste0("Jaccard Index vs Signal Strength | ", slice),
             subtitle = "Consistency of significant sets across rotations at each alpha level") +
        theme_bw() +
        theme(legend.position = "right")
      save_figure(p_jaccard_alpha, file.path(module2_fig_dir, "jaccard_by_alpha"),
                  width = 12, height = 5)

      cat("Module 2 outputs saved to", module2_dir, "\n\n")
    }

    # ==============================================================================
    # RUNTIME  (both modes - no alpha ground truth needed)
    # ==============================================================================

    cat("===== RUNTIME =====\n\n")

    runtime_summary <- all_data %>%
      select(tool, tool_label, angle, runtime_sec) %>%
      distinct() %>%
      arrange(tool, angle)
    write.csv(runtime_summary, file.path(module2_data_dir, "runtime.csv"), row.names = FALSE)

    p_runtime <- ggplot(runtime_summary,
                        aes(x = factor(angle), y = runtime_sec, fill = tool_label)) +
      geom_col(position = "dodge", color = "black", linewidth = 0.2) +
      scale_fill_manual(values = method_colors, name = "Method") +
      labs(x = "Rotation angle (\u00b0)", y = "Runtime (seconds)",
           title = paste0("Runtime per Angle - All Methods | ", slice, "/", mode)) +
      theme_bw() +
      theme(legend.position = "right")
    save_figure(p_runtime, file.path(module2_fig_dir, "runtime_barplot"),
                width = 10, height = 6)

    cat("Runtime outputs saved to", module2_dir, "\n")

    cat("\nAll metrics saved for ", slice, " / ", mode, ".\n")
  }
}

cat("\nAll slices/modes processed.\n")