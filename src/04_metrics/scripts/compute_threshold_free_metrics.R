#!/usr/bin/env Rscript
# ==============================================================================
# compute_threshold_free_metrics.R
# Threshold-free metrics for the SVG rotation-invariance benchmark.
#
# Design principle: the two modes are COMPLEMENTARY (no metric is computed in
# both modes, avoiding redundancy and biased set-overlap on the non-independent
# simulated features):
#
#   simulated -> auPRC + Kendall tau (per-angle, all 4)   [performance + invariance]
#   whole     -> top-K Jaccard + 4-way Venn (K = 2000)    [prioritization stability]
#
# Ranking scores (guide section 4 mapping; higher = more spatial signal):
#   moransi   -> I                  (identity)
#   spatialde -> FSV                (identity)
#   nnsvg     -> prop_sv            (identity)
#   sparkx    -> -log10(adjustedPval)
#   smash     -> -log10(pval)        (raw p-value)
#
# Metrics computed inline (no wrapper helpers): PRROC::pr.curve() for auPRC,
# cor(..., method="kendall") for Kendall tau.
#
# Input:  src/03_benchmark/outputs/{slice}/{mode}/{method}/scdesign3_angle{angle}_results.rds
# Output: src/04_metrics/outputs/{slice}/{mode}/{module2_threshold_free|module1_topk|runtime}/...
#         src/04_metrics/outputs/cross_slice/{simulated|whole}/figures/...
# ==============================================================================

# ==============================================================================
# SECTION 0: SETUP AND CONFIGURATION
# ==============================================================================

library(ggplot2)
library(reshape2)
library(dplyr)
library(ggvenn)
library(PRROC)

# --- Select which slices and modes to run ---
SLICES <- c("anterior1", "anterior2", "posterior1", "posterior2")
MODES  <- c("simulated", "whole")

# --- Top-K for whole-mode set overlap (guide section 4, Li et al. 2025) ---
TOP_K <- 2000

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# --- Score column and transform per method ---
# transform = "identity"  : use the column as-is (higher = more signal)
# transform = "neglog10_*": -log10(p-value) (higher = more signal)
score_config <- list(
  moransi   = list(col = "I",             transform = "identity"),
  spatialde = list(col = "FSV",           transform = "identity"),
  nnsvg     = list(col = "LR_stat",       transform = "identity"),
  sparkx    = list(col = "adjustedPval",  transform = "neglog10"),
  smash     = list(col = "pval",          transform = "neglog10")
)

# --- Method colors (Wong palette, colorblind-safe) ---
method_colors <- c(
  SPARKX    = "#0072B2",
  NNSVG     = "#D55E00",
  SPATIALDE = "#009E73",
  MORANSI   = "#CC79A7",
  SMASH     = "#F0E442"
)

angles  <- c(0, 30, 45, 60)
pairs_all <- list(c("0", "30"), c("0", "45"), c("0", "60"),
                  c("30", "45"), c("30", "60"), c("45", "60"))

# --- Angle colors for PR curves and scatter plots ---
angle_colors <- c(
  "0"  = "#000000",
  "30" = "#D55E00",
  "45" = "#0072B2",
  "60" = "#009E73"
)

# --- Helper: save both PNG (cairo) and PDF (reused from compute_all_metrics.R) ---
save_figure <- function(plot, filename_base, width, height, dpi = 300) {
  ggsave(paste0(filename_base, ".png"), plot, width = width, height = height,
         dpi = dpi, bg = "white")
  ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height, bg = "white")
}

# --- Helper: Jaccard index for two sets (reused from compute_all_metrics.R) ---
jaccard_index <- function(set_a, set_b) {
  inter <- length(intersect(set_a, set_b))
  uni   <- length(union(set_a, set_b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

# --- Containers for cross-slice aggregation ---
all_threshold_free <- list()
all_topk_jaccard   <- list()

# --- Container for PR curve data (simulated mode, per-slice) ---
all_pr_curves <- list()

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

    module_dir <- if (mode == "simulated") {
      file.path(metrics_root, "module2_threshold_free")
    } else {
      file.path(metrics_root, "module1_topk")
    }
    module_data_dir <- file.path(module_dir, "data")
    module_fig_dir  <- file.path(module_dir, "figures")
    runtime_dir     <- file.path(metrics_root, "runtime")
    runtime_data_dir <- file.path(runtime_dir, "data")
    runtime_fig_dir  <- file.path(runtime_dir, "figures")

    for (d in c(module_data_dir, module_fig_dir,
                runtime_data_dir, runtime_fig_dir)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }

    message("\n========== Threshold-free metrics: ", slice, " / ", mode, " ==========\n")

    # --------------------------------------------------------------------------
    # SECTION 1: DATA LOADING
    # --------------------------------------------------------------------------
    # Reads res_mtest from each method x angle, applies the ranking-score
    # transform inline (identity or -log10), and parses the alpha suffix from
    # feature names in simulated mode (gene_alpha; NA in whole mode).
    # --------------------------------------------------------------------------

    cat("Loading benchmark results for", slice, "/", mode, "...\n")

    tools <- list.files(benchmark_root)
    tools <- tools[tools %in% names(score_config)]

    all_data <- do.call(rbind, lapply(tools, function(tool) {
      cfg <- score_config[[tool]]
      if (is.null(cfg)) return(NULL)

      do.call(rbind, lapply(angles, function(angle) {
        rds_file <- file.path(benchmark_root, tool,
                              paste0("scdesign3_angle", angle, "_results.rds"))
        runtime_file <- file.path(benchmark_root, tool,
                                  paste0("scdesign3_angle", angle, "_runtime.csv"))
        if (!file.exists(rds_file)) return(NULL)

        df <- readRDS(rds_file)$res_mtest
        feature_names <- rownames(df)

        # Feature names: "gene_alpha" in simulated mode; plain gene in whole mode.
        parts <- strsplit(feature_names, "_")
        gene  <- sapply(parts, `[`, 1)
        alpha <- suppressWarnings(as.numeric(sapply(parts, function(x) {
          if (length(x) < 2) return(NA_real_)
          as.numeric(paste(x[-1], collapse = "_"))
        })))

        # Apply ranking-score transform inline (no helper function).
        # pmax with .Machine$double.xmin guards log10(0).
        if (cfg$transform == "identity") {
          score <- df[[cfg$col]]
        } else if (cfg$transform == "neglog10") {
          p <- pmax(df[[cfg$col]], .Machine$double.xmin)
          score <- -log10(p)
        } else {
          stop("Unknown transform '", cfg$transform, "' for tool '", tool, "'")
        }

        runtime_val <- if (file.exists(runtime_file)) {
          read.csv(runtime_file)$elapsed_sec
        } else {
          NA_real_
        }

        data.frame(
          tool = tool, angle = angle, gene = gene, alpha = alpha,
          feature = feature_names, score = score,
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
    } else {
      cat("Features per method x angle:", length(unique(all_data$feature)), "\n")
    }
    cat("\n")

    # ==============================================================================
    # SIMULATED MODE: MODULE 2 - THRESHOLD-FREE (auPRC + Kendall tau, per-angle)
    # ==============================================================================
    # auPRC: PRROC::pr.curve() with positives = (alpha > 0); auc.integral.
    # Kendall tau: cor(score, alpha, method="kendall") -- ordinal recovery of
    # the continuous signal-strength gradient.
    # Computed at ALL 4 angles so the cross-slice faceted heatmap (guide
    # section 5) can show performance + invariance simultaneously.
    # ==============================================================================

    if (mode == "simulated") {

      cat("===== MODULE 2: THRESHOLD-FREE (auPRC + Kendall tau) =====\n\n")

      threshold_free <- do.call(rbind, lapply(tools, function(tool) {
        do.call(rbind, lapply(angles, function(angle) {
          sub <- all_data[all_data$tool == tool & all_data$angle == angle, ]
          if (nrow(sub) == 0) return(NULL)

          pos <- sub$alpha > 0
          score <- sub$score

          # auPRC via PRROC::pr.curve (inline, no wrapper).
          # curve = TRUE so the PR curve coordinates are collected for the
          # per-slice PR curve panel below.
          pr <- pr.curve(scores.class0 = score[pos],
                         scores.class1 = score[!pos],
                         curve = TRUE)
          auprc <- pr$auc.integral

          # Kendall tau vs continuous alpha (inline, no wrapper).
          kendall <- cor(score, sub$alpha, method = "kendall",
                         use = "complete.obs")

          data.frame(
            slice = slice, tool = tool, tool_label = toupper(tool),
            angle = angle, auPRC = auprc, kendall_tau = kendall,
            n_pos = sum(pos), n_neg = sum(!pos),
            stringsAsFactors = FALSE
          )
        }))
      }))

      if (!is.null(threshold_free) && nrow(threshold_free) > 0) {
        write.csv(threshold_free,
                  file.path(module_data_dir, "threshold_free.csv"),
                  row.names = FALSE)

        cat("auPRC + Kendall tau (all angles):\n")
        print(threshold_free[, c("tool_label", "angle", "auPRC", "kendall_tau")])
        cat("\n")

        all_threshold_free[[length(all_threshold_free) + 1]] <- threshold_free
      }

      # --- Collect PR curve coordinates for the per-slice panel ---
      pr_curve_list <- do.call(rbind, lapply(tools, function(tool) {
        do.call(rbind, lapply(angles, function(angle) {
          sub <- all_data[all_data$tool == tool & all_data$angle == angle, ]
          if (nrow(sub) == 0) return(NULL)
          pos <- sub$alpha > 0
          pr <- pr.curve(scores.class0 = sub$score[pos],
                         scores.class1 = sub$score[!pos],
                         curve = TRUE)
          if (is.null(pr$curve)) return(NULL)
          data.frame(
            tool = tool, tool_label = toupper(tool),
            angle = as.character(angle),
            recall = pr$curve[, 1],
            precision = pr$curve[, 2],
            stringsAsFactors = FALSE
          )
        }))
      }))

      if (!is.null(pr_curve_list) && nrow(pr_curve_list) > 0) {
        pr_curve_list$angle <- factor(pr_curve_list$angle, levels = c("0","30","45","60"))

        p_prcurves <- ggplot(pr_curve_list,
                             aes(x = recall, y = precision, color = angle)) +
          geom_path(linewidth = 0.8) +
          facet_wrap(~ tool_label, ncol = 5) +
          scale_color_manual(values = angle_colors, name = "Angle (\u00b0)") +
          scale_x_continuous(limits = c(0, 1)) +
          scale_y_continuous(limits = c(0, 1)) +
          labs(x = "Recall", y = "Precision",
               title = paste0("Precision-Recall Curves Across Rotations | ", slice),
               subtitle = "Positives = alpha > 0 | Overlapping curves = rotation invariance") +
          theme_bw() +
          theme(strip.text = element_text(face = "bold", size = 11),
                legend.position = "bottom")
        save_figure(p_prcurves, file.path(module_fig_dir, "prcurves_panel"),
                    width = 16, height = 5)
      }

      # --- Kendall scatter panel: rank(score) vs alpha, faceted by method ---
      # Ranks computed within each tool x angle group; loess trend per angle.
      scatter_data <- do.call(rbind, lapply(tools, function(tool) {
        do.call(rbind, lapply(angles, function(angle) {
          sub <- all_data[all_data$tool == tool & all_data$angle == angle, ]
          if (nrow(sub) == 0) return(NULL)
          data.frame(
            tool = tool, tool_label = toupper(tool),
            angle = as.character(angle),
            alpha = sub$alpha,
            rank_score = rank(sub$score, ties.method = "average"),
            stringsAsFactors = FALSE
          )
        }))
      }))

      if (!is.null(scatter_data) && nrow(scatter_data) > 0) {
        scatter_data$angle <- factor(scatter_data$angle, levels = c("0","30","45","60"))

        p_scatter <- ggplot(scatter_data,
                            aes(x = alpha, y = rank_score, color = angle)) +
          geom_point(alpha = 0.3, size = 0.8) +
          geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
          facet_wrap(~ tool_label, ncol = 5) +
          scale_color_manual(values = angle_colors, name = "Angle (\u00b0)") +
          labs(x = "Signal strength (alpha)", y = "Rank(score)",
               title = paste0("Kendall \u03c4 Scatter: Rank vs Signal Gradient | ", slice),
               subtitle = "Monotonic step-up = good ordinal recovery | Angle overlap = invariance") +
          theme_bw() +
          theme(strip.text = element_text(face = "bold", size = 11),
                legend.position = "bottom")
        save_figure(p_scatter, file.path(module_fig_dir, "kendall_scatter_panel"),
                    width = 16, height = 5)
      }

      cat("Module 2 (threshold-free) outputs saved to", module_dir, "\n\n")
    }

    # ==============================================================================
    # WHOLE MODE: MODULE 1 - TOP-K SET OVERLAP (Jaccard + Venn, K = 2000)
    # ==============================================================================
    # Genes ranked by the continuous ranking score (higher = more spatial); the
    # top-K (K = 2000) features per method x angle define the sets. Jaccard is
    # computed over all 6 angle pairs; Venn diagrams are 4-way (0/30/45/60).
    # Homogeneous set sizes isolate prioritization stability under rotation
    # (guide section 4, Li et al. 2025).
    # ==============================================================================

    if (mode == "whole") {

      cat("===== MODULE 1: TOP-K SET OVERLAP (K = ", TOP_K, ") =====\n\n")

      topk_sets_by_tool <- lapply(tools, function(tool) {
        sets <- lapply(angles, function(angle) {
          sub <- all_data[all_data$tool == tool & all_data$angle == angle, ]
          if (nrow(sub) == 0) return(character(0))
          # Rank by score descending; take top-K features.
          sub <- sub[order(sub$score, decreasing = TRUE), ]
          head(sub$feature, TOP_K)
        })
        names(sets) <- as.character(angles)
        sets
      })
      names(topk_sets_by_tool) <- tools

      jaccard_topk <- do.call(rbind, lapply(tools, function(tool) {
        sets <- topk_sets_by_tool[[tool]]
        do.call(rbind, lapply(pairs_all, function(pair) {
          a <- pair[1]; b <- pair[2]
          if (!a %in% names(sets) || !b %in% names(sets)) return(NULL)
          data.frame(
            slice = slice, tool = tool, tool_label = toupper(tool),
            angle_a = as.numeric(a), angle_b = as.numeric(b),
            jaccard = jaccard_index(sets[[a]], sets[[b]]),
            pair_label = paste0(a, "\u00b0 vs ", b, "\u00b0"),
            K = TOP_K,
            stringsAsFactors = FALSE
          )
        }))
      }))

      if (!is.null(jaccard_topk) && nrow(jaccard_topk) > 0) {
        write.csv(jaccard_topk, file.path(module_data_dir, "jaccard_topk.csv"),
                  row.names = FALSE)

        cat("Top-K (K=", TOP_K, ") Jaccard across angle pairs:\n")
        print(jaccard_topk[, c("tool_label", "pair_label", "jaccard")])
        cat("\n")

        all_topk_jaccard[[length(all_topk_jaccard) + 1]] <- jaccard_topk
      }

      # --- Per-slice Jaccard heatmap (methods x angle pairs) ---
      if (!is.null(jaccard_topk) && nrow(jaccard_topk) > 0) {
        p_jaccard <- ggplot(jaccard_topk,
                            aes(x = pair_label, y = tool_label, fill = jaccard)) +
          geom_tile(color = "white", linewidth = 0.8) +
          geom_text(aes(label = sprintf("%.3f", jaccard)),
                    size = 3.5, fontface = "bold") +
          scale_fill_gradient2(low = "#D73027", mid = "#FFFFCC", high = "#1A9850",
                               midpoint = 0.95, limits = c(0, 1),
                               name = "Jaccard index") +
          labs(x = "Angle pair", y = "Method",
               title = paste0("Top-", TOP_K, " Set Overlap: Jaccard Across Angle Pairs | ", slice),
               subtitle = paste0("K = ", TOP_K,
                                 " top-ranked genes per method x angle")) +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 30, hjust = 1))
        save_figure(p_jaccard, file.path(module_fig_dir, "jaccard_topk_heatmap"),
                    width = 10, height = 6)
      }

      # --- 4-way Venn diagrams on top-K sets (one per method) ---
      for (tool in tools) {
        sets <- topk_sets_by_tool[[tool]]
        if (length(sets) < 4) next
        venn_data <- list(
          Original     = sets[["0"]],
          `30 degree`  = sets[["30"]],
          `45 degree`  = sets[["45"]],
          `60 degree`  = sets[["60"]]
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
          labs(title = paste0(toupper(tool), ": Top-", TOP_K,
                              " SVG Overlap Across Rotations (", slice, ")"),
               subtitle = paste0("Top-", TOP_K,
                                 " by ranking score | Intersection = stable prioritization")) +
          theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
        save_figure(p_venn, file.path(module_fig_dir, paste0("venn_topk_", tool)),
                    width = 10, height = 8)
      }

      cat("Module 1 (top-K) outputs saved to", module_dir, "\n\n")
    }

    # ==============================================================================
    # RUNTIME (both modes)
    # ==============================================================================

    cat("===== RUNTIME =====\n\n")

    runtime_summary <- all_data %>%
      select(tool, tool_label, angle, runtime_sec) %>%
      distinct() %>%
      arrange(tool, angle)
    write.csv(runtime_summary, file.path(runtime_data_dir, "runtime.csv"),
              row.names = FALSE)

    p_runtime <- ggplot(runtime_summary,
                        aes(x = factor(angle), y = runtime_sec, fill = tool_label)) +
      geom_col(position = "dodge", color = "black", linewidth = 0.2) +
      scale_fill_manual(values = method_colors, name = "Method") +
      labs(x = "Rotation angle (\u00b0)", y = "Runtime (seconds)",
           title = paste0("Runtime per Angle - All Methods | ", slice, "/", mode)) +
      theme_bw() +
      theme(legend.position = "right")
    save_figure(p_runtime, file.path(runtime_fig_dir, "runtime_barplot"),
                width = 10, height = 6)

    cat("Runtime outputs saved to", runtime_dir, "\n")
    cat("\nAll metrics saved for ", slice, " / ", mode, ".\n")
  }
}

# ==============================================================================
# SECTION 5: CROSS-SLICE FLAGSHIP FIGURES (guide section 5)
# ==============================================================================
# Faceted heatmaps consolidating all slices: rows = methods, columns = angles
# (or angle pairs), faceted horizontally by slice. This is the guide's
# section-5 architecture showing performance + invariance in a single view.
# ==============================================================================

cat("\n===== CROSS-SLICE FLAGSHIP FIGURES =====\n\n")

cross_fig_root <- file.path(project_root, "src", "04_metrics", "outputs", "cross_slice")
cross_sim_dir  <- file.path(cross_fig_root, "simulated", "figures")
cross_whole_dir <- file.path(cross_fig_root, "whole", "figures")
dir.create(cross_sim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cross_whole_dir, recursive = TRUE, showWarnings = FALSE)

# --- Simulated: auPRC + Kendall tau faceted heatmaps ---

if (length(all_threshold_free) > 0) {
  tf_all <- do.call(rbind, all_threshold_free)
  tf_all$slice <- factor(tf_all$slice, levels = SLICES)

  # auPRC faceted heatmap (methods x angles, faceted by slice)
  p_auprc_heat <- ggplot(tf_all,
                         aes(x = factor(angle), y = tool_label, fill = auPRC)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", auPRC)), size = 3, fontface = "bold") +
    scale_fill_gradient2(low = "#D73027", mid = "#FFFFBF", high = "#1A9850",
                         midpoint = 0.95, limits = c(0.9, 1), name = "auPRC") +
    facet_wrap(~ slice, nrow = 1) +
    labs(x = "Rotation angle (\u00b0)", y = "Method",
         title = "auPRC Across Rotations and Slices (Simulated)",
         subtitle = "Stable color block = rotation invariance | Fading at oblique angles = vulnerability") +
    theme_bw() +
    theme(strip.text = element_text(face = "bold", size = 11),
          axis.text.x = element_text(size = 10))
  save_figure(p_auprc_heat, file.path(cross_sim_dir, "auprc_faceted_heatmap"),
              width = 16, height = 5)

  # Kendall tau faceted heatmap (diverging palette centered at 0)
  p_kendall_heat <- ggplot(tf_all,
                           aes(x = factor(angle), y = tool_label, fill = kendall_tau)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", kendall_tau)), size = 3, fontface = "bold") +
    scale_fill_gradient2(low = "#762A83", mid = "#F7F7F7", high = "#1B7837",
                         midpoint = 0, limits = c(-1, 1), name = "Kendall \u03c4") +
    facet_wrap(~ slice, nrow = 1) +
    labs(x = "Rotation angle (\u00b0)", y = "Method",
         title = "Kendall \u03c4 vs Signal Gradient Across Rotations and Slices (Simulated)",
         subtitle = "Ordinal recovery of alpha | Stable block = invariance") +
    theme_bw() +
    theme(strip.text = element_text(face = "bold", size = 11),
          axis.text.x = element_text(size = 10))
  save_figure(p_kendall_heat, file.path(cross_sim_dir, "kendall_faceted_heatmap"),
              width = 16, height = 5)

  cat("Simulated cross-slice flagship heatmaps saved to", cross_sim_dir, "\n")
} else {
  cat("No simulated threshold-free data collected -- skipping simulated cross-slice figures.\n")
}

# --- Whole: top-K Jaccard faceted heatmap ---

if (length(all_topk_jaccard) > 0) {
  j_all <- do.call(rbind, all_topk_jaccard)
  j_all$slice <- factor(j_all$slice, levels = SLICES)

  p_jaccard_heat <- ggplot(j_all,
                           aes(x = pair_label, y = tool_label, fill = jaccard)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", jaccard)), size = 3, fontface = "bold") +
    scale_fill_gradient2(low = "#D73027", mid = "#FFFFCC", high = "#1A9850",
                         midpoint = 0.95, limits = c(0, 1), name = "Jaccard") +
    facet_wrap(~ slice, nrow = 1) +
    labs(x = "Angle pair", y = "Method",
         title = paste0("Top-", TOP_K, " Jaccard Across Rotations and Slices (Whole)"),
         subtitle = paste0("Prioritization stability under rotation | K = ", TOP_K,
                           " top-ranked genes")) +
    theme_bw() +
    theme(strip.text = element_text(face = "bold", size = 11),
          axis.text.x = element_text(angle = 30, hjust = 1))
  save_figure(p_jaccard_heat, file.path(cross_whole_dir, "jaccard_topk_faceted_heatmap"),
              width = 16, height = 5)

  cat("Whole cross-slice flagship heatmap saved to", cross_whole_dir, "\n")
} else {
  cat("No whole top-K Jaccard data collected -- skipping whole cross-slice figure.\n")
}

cat("\nAll slices/modes processed.\n")
