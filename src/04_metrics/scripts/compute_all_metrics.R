#!/usr/bin/env Rscript
# ==============================================================================
# compute_all_metrics.R
# Comprehensive metrics computation for SVG rotation invariance benchmark.
#
# Strictly separated into:
#   MODULE 1: Evaluating Geometric Robustness (Rotation Invariance)  — all angles
#   MODULE 2: Evaluating Statistical Performance (Simulation)         — 0° only
#   SUPPLEMENTARY: Extra metrics (AUC-ROC, classification, confusion)  — kept for reference
# ==============================================================================

# ==============================================================================
# SECTION 0: SETUP AND CONFIGURATION
# ==============================================================================

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(tidyr)

ggVennDiagram_available <- requireNamespace("ggVennDiagram", quietly = TRUE)
if (ggVennDiagram_available) {
  library(ggVennDiagram)
}

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

benchmark_root <- file.path(project_root, "src", "03_benchmark", "outputs")
metrics_root   <- file.path(project_root, "src", "04_metrics", "outputs")

module1_dir       <- file.path(metrics_root, "module1_rotation")
module1_data_dir  <- file.path(module1_dir, "data")
module1_fig_dir   <- file.path(module1_dir, "figures")
module2_dir       <- file.path(metrics_root, "module2_performance")
module2_data_dir  <- file.path(module2_dir, "data")
module2_fig_dir   <- file.path(module2_dir, "figures")
supp_dir          <- file.path(metrics_root, "supplementary")
supp_data_dir     <- file.path(supp_dir, "data")
supp_fig_dir      <- file.path(supp_dir, "figures")

for (d in c(module1_data_dir, module1_fig_dir,
            module2_data_dir, module2_fig_dir,
            supp_data_dir, supp_fig_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

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
  moransi   = "pval_norm_fdr_bh",
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

# --- Bioinformatics publication-ready theme ---
theme_bioinformatics <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major      = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor      = element_blank(),
      panel.border          = element_rect(color = "black", linewidth = 0.5),
      strip.background      = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
      strip.text            = element_text(face = "bold", size = base_size),
      axis.title            = element_text(face = "bold", size = base_size + 1),
      axis.text             = element_text(size = base_size),
      legend.title          = element_text(face = "bold", size = base_size),
      legend.text           = element_text(size = base_size - 1),
      legend.background     = element_rect(fill = "white", color = "gray80"),
      plot.title            = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      plot.subtitle         = element_text(size = base_size, hjust = 0.5, color = "gray30"),
      plot.caption          = element_text(size = base_size - 2, hjust = 1, color = "gray50"),
      plot.margin           = margin(10, 10, 10, 10)
    )
}

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

# --- Helper: AUPRC via trapezoidal integration ---
compute_auprc <- function(truth, score) {
  truth <- as.integer(truth); score <- as.numeric(score)
  if (length(truth) == 0L || all(is.na(truth)) || length(unique(truth)) < 2L) return(NA_real_)
  valid_idx <- which(!is.na(truth) & !is.na(score))
  truth <- truth[valid_idx]; score <- score[valid_idx]
  if (!any(truth == 1L) || !any(truth == 0L)) return(NA_real_)
  ordering <- order(score, decreasing = TRUE); truth <- truth[ordering]
  cumulative_positives <- cumsum(truth == 1L)
  precision <- cumulative_positives / seq_along(truth)
  recall <- cumulative_positives / sum(truth == 1L)
  recall_previous <- c(0, recall[-length(recall)])
  sum((recall - recall_previous) * precision)
}

# --- Helper: AU-ROC ---
compute_auroc <- function(truth, score) {
  truth <- as.integer(truth); score <- as.numeric(score)
  if (length(truth) == 0L || all(is.na(truth)) || length(unique(truth)) < 2L) return(NA_real_)
  valid_idx <- which(!is.na(truth) & !is.na(score))
  truth <- truth[valid_idx]; score <- score[valid_idx]
  if (!any(truth == 1L) || !any(truth == 0L)) return(NA_real_)
  n_pos <- sum(truth == 1L); n_neg <- sum(truth == 0L)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[truth == 1L]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

# --- Helper: classification metrics ---
compute_classification_metrics <- function(truth, predicted) {
  tp <- as.numeric(sum(truth == 1 & predicted == 1, na.rm = TRUE))
  tn <- as.numeric(sum(truth == 0 & predicted == 0, na.rm = TRUE))
  fp <- as.numeric(sum(truth == 0 & predicted == 1, na.rm = TRUE))
  fn <- as.numeric(sum(truth == 1 & predicted == 0, na.rm = TRUE))
  sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
  specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_)
  precision   <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
  fdr         <- ifelse((tp + fp) > 0, fp / (tp + fp), NA_real_)
  f1          <- ifelse((precision + sensitivity) > 0,
                        2 * precision * sensitivity / (precision + sensitivity), NA_real_)
  mcc_num <- (tp * tn - fp * fn)
  mcc_den <- sqrt(max(0, (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
  mcc <- ifelse(mcc_den > 0, mcc_num / mcc_den, NA_real_)
  balanced_acc <- (sensitivity + specificity) / 2
  data.frame(TP = tp, TN = tn, FP = fp, FN = fn,
             Sensitivity = sensitivity, Specificity = specificity,
             Precision = precision, FDR = fdr, F1 = f1,
             MCC = mcc, BalancedAccuracy = balanced_acc,
             stringsAsFactors = FALSE)
}

# ==============================================================================
# SECTION 1: DATA LOADING — unified load into a single data frame
# ==============================================================================

cat("Loading benchmark results...\n")

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

    parts <- strsplit(feature_names, "_")
    gene  <- sapply(parts, `[`, 1)
    alpha <- as.numeric(sapply(parts, function(x) paste(x[-1], collapse = "_")))

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

all_data$tool_label <- toupper(all_data$tool)
all_data$angle <- as.numeric(all_data$angle)

cat("Loaded", nrow(all_data), "rows for",
    length(unique(all_data$tool)), "methods across",
    length(unique(all_data$angle)), "angles.\n")
cat("Methods:", paste(unique(all_data$tool_label), collapse = ", "), "\n")
cat("Alpha levels:", paste(sort(unique(all_data$alpha)), collapse = ", "), "\n\n")

# --- Classification metrics at adj p < 0.05 (all angles, needed by M2 + supp) ---
class_df <- do.call(rbind, lapply(tools, function(tool) {
  do.call(rbind, lapply(angles, function(angle) {
    sub <- all_data[all_data$tool == tool & all_data$angle == angle, ]
    if (nrow(sub) == 0) return(NULL)
    truth <- as.integer(sub$alpha > 0)
    predicted <- as.integer(sub$significant)
    cm <- compute_classification_metrics(truth, predicted)
    cm$tool <- tool; cm$tool_label <- toupper(tool); cm$angle <- angle
    cm
  }))
}))
write.csv(class_df, file.path(supp_data_dir, "classification_metrics.csv"), row.names = FALSE)


# ==============================================================================
### MODULE 1: ROTATION INVARIANCE ###
# ==============================================================================
# Goal: Evaluate if methods yield the same results when spatial coordinates
# are rotated (0, 30, 45, 60 degrees).
# Data subsetting: ALL angles.
# ==============================================================================

cat("===== MODULE 1: ROTATION INVARIANCE =====\n\n")

# --- Metric 1.1: Set Overlap (Jaccard & Venn) -------------------------------

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
write.csv(jaccard_results, file.path(module1_data_dir, "jaccard.csv"), row.names = FALSE)

# Venn diagrams (4-way: 0, 30, 45, 60) per method
# Distinct border colors per angle; low-alpha fill so intersection counts are readable
angle_edge_colors <- c("1" = "#0072B2", "2" = "#D55E00", "3" = "#009E73", "4" = "#CC79A7")
angle_labels <- c("1" = "0\u00b0", "2" = "30\u00b0", "3" = "45\u00b0", "4" = "60\u00b0")

if (ggVennDiagram_available) {
  for (tool in tools) {
    sets <- sig_sets_by_tool[[tool]]
    if (length(sets) < 4) next
    venn_data <- list(`0` = sets[["0"]], `30` = sets[["30"]],
                      `45` = sets[["45"]], `60` = sets[["60"]])
    p_venn <- ggVennDiagram(venn_data, label = "count", label_alpha = 0,
                            edge_size = 0.8)
    # Override the edge colour mapping to use distinct colors per angle
    p_venn$layers[[2]]$mapping <- aes(colour = factor(id), group = id,
                                       linetype = I(linetype), linewidth = I(linewidth))
    p_venn <- p_venn +
      scale_fill_gradient(low = "#F7F7F7", high = method_colors[toupper(tool)],
                          name = "Count", aesthetics = "fill") +
      scale_colour_manual(values = angle_edge_colors, name = "Angle",
                          labels = angle_labels) +
      labs(title = paste0(toupper(tool), ": Significant SVG Overlap Across Rotations"),
           subtitle = "Adjusted p < 0.05 | Intersection = consistent across all angles") +
      theme_bioinformatics() +
      theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
            legend.position = "right")
    save_figure(p_venn, file.path(module1_fig_dir, paste0("venn_", tool)),
                width = 10, height = 8)
  }
} else {
  message("ggVennDiagram not available - skipping Venn diagrams")
}

# Cross-method Jaccard heatmap
p_jaccard <- ggplot(jaccard_results,
                    aes(x = pair_label, y = tool_label, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", jaccard)), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#D73027", mid = "#FFFFCC", high = "#1A9850",
                       midpoint = 0.95, limits = c(0, 1), name = "Jaccard index") +
  labs(x = "Angle pair", y = "Method",
       title = "Set Overlap: Jaccard Index Across Angle Pairs",
       subtitle = "Jaccard = |intersection| / |union| of significant gene sets") +
  theme_bioinformatics() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_figure(p_jaccard, file.path(module1_fig_dir, "jaccard_heatmap"), width = 10, height = 6)

# --- Metric 1.2: Score/Statistic Consistency (Plot & Table) ----------------

score_wide <- all_data %>%
  select(tool, tool_label, feature, angle, score) %>%
  pivot_wider(
    id_cols = c(tool, tool_label, feature),
    names_from = angle,
    values_from = score,
    names_prefix = "angle_"
  )

correlation_summary <- score_wide %>%
  group_by(tool, tool_label) %>%
  summarise(
    pearson_30  = cor(angle_0, angle_30,  method = "pearson",  use = "complete.obs"),
    spearman_30 = cor(angle_0, angle_30,  method = "spearman", use = "complete.obs"),
    pearson_45  = cor(angle_0, angle_45,  method = "pearson",  use = "complete.obs"),
    spearman_45 = cor(angle_0, angle_45,  method = "spearman", use = "complete.obs"),
    pearson_60  = cor(angle_0, angle_60,  method = "pearson",  use = "complete.obs"),
    spearman_60 = cor(angle_0, angle_60,  method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )
write.csv(correlation_summary,
          file.path(module1_data_dir, "correlation_summary.csv"), row.names = FALSE)

cat("Correlation summary (0 vs rotated):\n")
print(correlation_summary)

# Faceted scatter plot per method: 0 deg vs each rotated angle
# Subtitle makes clear what statistic is being plotted (score column + type)
scatter_df <- do.call(rbind, lapply(tools, function(tool) {
  sub <- score_wide[score_wide$tool == tool, ]
  do.call(rbind, lapply(pairs_vs_0, function(pair) {
    a <- pair[1]; b <- pair[2]
    col_a <- paste0("angle_", a); col_b <- paste0("angle_", b)
    if (!col_a %in% names(sub) || !col_b %in% names(sub)) return(NULL)
    data.frame(
      tool = tool, tool_label = toupper(tool),
      score_0 = sub[[col_a]], score_theta = sub[[col_b]],
      angle_a = a, angle_b = b,
      stringsAsFactors = FALSE
    )
  }))
}))

for (tool in tools) {
  sub <- scatter_df[scatter_df$tool == tool, ]
  if (nrow(sub) == 0) next
  cfg <- score_config[[tool]]
  stat_label <- if (cfg$negate)
    paste0(cfg$col, " (negated p-value)")
  else
    paste0(cfg$col, " (test statistic)")

  sub$label <- factor(paste0(sub$angle_a, "\u00b0 vs ", sub$angle_b, "\u00b0"),
                      levels = paste0("0\u00b0 vs ", c("30\u00b0", "45\u00b0", "60\u00b0")))
  lims <- range(c(sub$score_0, sub$score_theta), na.rm = TRUE)
  p <- ggplot(sub, aes(x = score_0, y = score_theta)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(alpha = 0.4, size = 1.2, color = method_colors[toupper(tool)]) +
    facet_wrap(~ label, ncol = 3) +
    coord_fixed(xlim = lims, ylim = lims) +
    labs(x = "Statistic at 0\u00b0", y = "Statistic at theta\u00b0",
         title = paste0(toupper(tool), ": Score Consistency Across Rotations"),
         subtitle = paste0("Score: ", stat_label, " | Diagonal = perfect rotation invariance")) +
    theme_bioinformatics()
  save_figure(p, file.path(module1_fig_dir, paste0("scatter_facet_", tool)),
              width = 12, height = 5)
}

# Rotation tau (computed here, saved to supplementary, used by supp trade-off plot)
rotation_tau <- do.call(rbind, lapply(tools, function(tool) {
  sub <- score_wide[score_wide$tool == tool, ]
  do.call(rbind, lapply(pairs_all, function(pair) {
    a <- pair[1]; b <- pair[2]
    col_a <- paste0("angle_", a); col_b <- paste0("angle_", b)
    if (!col_a %in% names(sub) || !col_b %in% names(sub)) return(NULL)
    tau <- cor(sub[[col_a]], sub[[col_b]], method = "kendall", use = "complete.obs")
    data.frame(tool = tool, tool_label = toupper(tool),
               angle_a = as.numeric(a), angle_b = as.numeric(b),
               tau_rotation = tau,
               pair_label = paste0(a, "\u00b0 vs ", b, "\u00b0"),
               stringsAsFactors = FALSE)
  }))
}))
write.csv(rotation_tau, file.path(supp_data_dir, "rotation_tau.csv"), row.names = FALSE)

rotation_summary <- rotation_tau %>%
  group_by(tool_label) %>%
  summarise(mean_tau = mean(tau_rotation, na.rm = TRUE),
            min_tau  = min(tau_rotation,  na.rm = TRUE),
            max_tau  = max(tau_rotation,  na.rm = TRUE)) %>%
  arrange(desc(mean_tau))
write.csv(rotation_summary,
          file.path(supp_data_dir, "rotation_invariance_summary.csv"), row.names = FALSE)

cat("Module 1 outputs saved to", module1_dir, "\n\n")


# ==============================================================================
### MODULE 2: STATISTICAL PERFORMANCE ###
# ==============================================================================
# Goal: Evaluate statistical power, calibration, and accuracy under simulated
# ground truths.
# Data subsetting: ONLY the 0° baseline rotation. No averaging across rotations.
# Ground truth: alpha = 0 is pure noise (TN); alpha > 0 is spatial signal (TP).
# ==============================================================================

cat("===== MODULE 2: STATISTICAL PERFORMANCE (0° only) =====\n\n")

baseline_data <- all_data %>% filter(angle == 0)

# --- Metric 2.1: Statistical Calibration (FPR) ------------------------------

# Isolate alpha == 0 (true negatives). FPR = fraction called significant.
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
       title = "Statistical Calibration: FPR at alpha = 0 (0\u00b0 baseline)",
       subtitle = "Dashed line = nominal 0.05 threshold | FPR >> 0.05 = overly aggressive") +
  theme_bioinformatics()
save_figure(p_fpr, file.path(module2_fig_dir, "fpr_barplot"), width = 9, height = 6)

# --- Metric 2.2: Limit of Detection (Sensitivity across signal gradients) ---

# Isolate alpha > 0 (true positives). Group by specific alpha level (dynamic).
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
       title = "Limit of Detection: Sensitivity Across Signal Gradients (0\u00b0)",
       subtitle = "At what signal strength does each method recognize an SVG?") +
  theme_bioinformatics() +
  theme(legend.position = "right")
save_figure(p_sens, file.path(module2_fig_dir, "sensitivity_lines"), width = 10, height = 6)

# --- Metric 2.3: Ranking Accuracy (Kendall's tau) ---------------------------

# Kendall tau between continuous alpha (0 to 1) and predicted score, at 0° only.
kendall_tau <- baseline_data %>%
  group_by(tool, tool_label) %>%
  summarise(
    kendall_tau = cor(alpha, score, method = "kendall", use = "complete.obs"),
    .groups = "drop"
  ) %>%
  arrange(desc(kendall_tau))
write.csv(kendall_tau, file.path(module2_data_dir, "kendall_tau.csv"), row.names = FALSE)

cat("\nKendall's tau (alpha vs score, 0° only):\n")
print(kendall_tau)

p_kendall <- ggplot(kendall_tau,
                    aes(x = reorder(tool_label, kendall_tau), y = kendall_tau,
                        fill = tool_label)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", kendall_tau)), vjust = -0.5, size = 4,
            fontface = "bold") +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(x = "Method", y = "Kendall's tau",
       title = "Ranking Accuracy: Kendall's tau at 0\u00b0",
       subtitle = "Correlation between true alpha and predicted spatial score") +
  theme_bioinformatics()
save_figure(p_kendall, file.path(module2_fig_dir, "kendall_barplot"), width = 9, height = 6)

# --- Confusion Matrix at 0 deg (moved from Supplementary) ------------------

conf_0 <- class_df[class_df$angle == 0, ]
conf_melt <- melt(conf_0, id.vars = c("tool_label", "angle"),
                  measure.vars = c("TP", "TN", "FP", "FN"),
                  variable.name = "Cell", value.name = "Count")

p_conf <- ggplot(conf_melt, aes(x = Cell, y = reorder(tool_label, Count), fill = Count)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Count), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "#FEE8C8", high = "#E6550D", name = "Count") +
  labs(x = "Confusion Matrix Cell", y = "Method",
       title = "Confusion Matrix at 0\u00b0 - All Methods",
       subtitle = "Adjusted p < 0.05 significance threshold") +
  theme_bioinformatics() +
  theme(axis.text.y = element_text(face = "bold", size = 11))
save_figure(p_conf, file.path(module2_fig_dir, "confusion_matrix_heatmap"),
            width = 10, height = 6)

# --- Runtime comparison (kept on main metrics, all angles) ------------------

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
       title = "Runtime per Angle - All Methods") +
  theme_bioinformatics() +
  theme(legend.position = "right")
save_figure(p_runtime, file.path(module2_fig_dir, "runtime_barplot"), width = 10, height = 6)

cat("Module 2 outputs saved to", module2_dir, "\n\n")


# ==============================================================================
### SUPPLEMENTARY METRICS ###
# ==============================================================================
# Extra metrics kept for reference: AUC-ROC, classification metrics by angle,
# rank stability, rotation invariance heatmaps/summary, detection-vs-invariance
# trade-off.
# ==============================================================================

cat("===== SUPPLEMENTARY METRICS =====\n\n")

# AUC-ROC / AUPRC per method per angle
supp_auroc <- all_data %>%
  group_by(tool, tool_label, angle) %>%
  summarise(
    auroc_any = compute_auroc(as.integer(alpha > 0), score),
    auprc_any = compute_auprc(as.integer(alpha > 0), score),
    .groups = "drop"
  )
write.csv(supp_auroc, file.path(supp_data_dir, "auroc_auprc.csv"), row.names = FALSE)

p_all_auroc <- ggplot(supp_auroc, aes(x = angle, y = auroc_any,
                                      color = tool_label, group = tool)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3) +
  scale_color_manual(values = method_colors, name = "Method") +
  labs(x = "Rotation angle (\u00b0)", y = "AUC-ROC",
       title = "AUC-ROC vs Rotation Angle - All Methods",
       subtitle = "Detection of spatially variable genes (alpha > 0, baseline = 0.5)") +
  theme_bioinformatics() +
  theme(legend.position = "right")
save_figure(p_all_auroc, file.path(supp_fig_dir, "auroc_by_method"), width = 10, height = 6)

# Classification metrics by angle (line plot, faceted by metric)
class_melt <- melt(class_df,
                   id.vars = c("tool_label", "angle"),
                   measure.vars = c("Sensitivity", "Specificity", "Precision",
                                    "FDR", "F1", "MCC", "BalancedAccuracy"),
                   variable.name = "Metric", value.name = "Value")
class_melt$angle <- factor(class_melt$angle, levels = as.character(angles))

p_class_angle <- ggplot(class_melt, aes(x = angle, y = Value,
                                        color = tool_label, group = tool_label)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  facet_wrap(~ Metric, ncol = 4, scales = "free_y") +
  scale_color_manual(values = method_colors, name = "Method") +
  labs(x = "Rotation angle (\u00b0)", y = "Metric value",
       title = "Classification Metrics vs Rotation Angle - All Methods",
       subtitle = "Adjusted p < 0.05 | Stable lines = rotation-invariant classification") +
  theme_bioinformatics() +
  theme(legend.position = "bottom")
save_figure(p_class_angle, file.path(supp_fig_dir, "classification_metrics_by_angle"),
            width = 14, height = 10)

# Rank stability facet (moved from Module 1)
rank_0_60 <- all_data %>%
  filter(angle %in% c(0, 60)) %>%
  select(tool, tool_label, feature, angle, score) %>%
  pivot_wider(id_cols = c(tool, tool_label, feature),
              names_from = angle, values_from = score,
              names_prefix = "angle_") %>%
  mutate(rank_0 = rank(-angle_0, ties.method = "average"),
         rank_60 = rank(-angle_60, ties.method = "average"))

p_rank <- ggplot(rank_0_60, aes(x = rank_0, y = rank_60)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.3, size = 0.6, aes(color = tool_label)) +
  facet_wrap(~ tool_label, ncol = 3) +
  scale_color_manual(values = method_colors, guide = "none") +
  coord_fixed() +
  labs(x = "Rank at 0\u00b0", y = "Rank at 60\u00b0",
       title = "Rank Stability: 0\u00b0 vs 60\u00b0 - All Methods",
       subtitle = "Diagonal = perfect rotation invariance; scatter = rotation-variant") +
  theme_bioinformatics()
save_figure(p_rank, file.path(supp_fig_dir, "rank_stability_facet"), width = 12, height = 9)

# Rotation invariance heatmap (moved from Module 1)
p_heat <- ggplot(rotation_tau, aes(x = pair_label, y = tool_label, fill = tau_rotation)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", tau_rotation)), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#D73027", mid = "#FFFFCC", high = "#1A9850",
                       midpoint = 0.9, limits = c(0.5, 1), name = "tau (rotation)") +
  labs(x = "Angle pair", y = "Method",
       title = "Rotation Invariance: Score Consistency Across Angle Pairs",
       subtitle = "Higher tau = more rotation-invariant") +
  theme_bioinformatics() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_figure(p_heat, file.path(supp_fig_dir, "rotation_invariance_heatmap"),
            width = 10, height = 6)

# Rotation invariance summary bar chart (moved from Module 1)
p_rot_summary <- ggplot(rotation_summary,
                        aes(x = reorder(tool_label, mean_tau), y = mean_tau)) +
  geom_col(fill = "#1A9850", alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = min_tau, ymax = max_tau), width = 0.2,
                color = "gray30", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.4f", mean_tau)), vjust = -0.5, size = 3.5,
            fontface = "bold") +
  labs(x = "Method", y = "Mean tau (rotation) (across all pairs)",
       title = "Rotation Invariance Summary",
       subtitle = "Higher = more consistent scores across rotations (error bars: min/max)") +
  coord_flip(ylim = c(0.5, 1)) +
  theme_bioinformatics()
save_figure(p_rot_summary, file.path(supp_fig_dir, "rotation_invariance_summary"),
            width = 9, height = 6)

# Detection vs Invariance trade-off
auroc_0 <- supp_auroc[supp_auroc$angle == 0, c("tool_label", "auroc_any")]
colnames(auroc_0)[2] <- "auroc_0"
tradeoff <- merge(auroc_0, rotation_summary, by = "tool_label")

p_tradeoff <- ggplot(tradeoff, aes(x = auroc_0, y = mean_tau, color = tool_label)) +
  geom_point(size = 5) +
  geom_text(aes(label = tool_label), vjust = -1.3, hjust = 0.5, size = 4,
            fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = method_colors, guide = "none") +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  labs(x = "AUC-ROC at 0\u00b0 (detection ability)",
       y = expression("Mean " * tau[rotation] ~ "(rotation invariance)"),
       title = "Detection vs Rotation Invariance Trade-off",
       subtitle = "Top-right = ideal: detects SVGs AND is rotation-invariant") +
  coord_cartesian(xlim = c(0.5, 1.01), ylim = c(0.5, 1.01)) +
  theme_bioinformatics()
save_figure(p_tradeoff, file.path(supp_fig_dir, "detection_vs_invariance"),
            width = 9, height = 8)

cat("Supplementary outputs saved to", supp_dir, "\n\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("===== SUMMARY =====\n")
cat("\nModule 1 - Rotation Invariance ranking (mean tau_rotation):\n")
print(rotation_summary)
cat("\nModule 2 - FPR at alpha=0 (0 deg):\n")
print(fpr_summary[, c("tool_label", "fpr")])
cat("\nModule 2 - Kendall's tau at 0 deg:\n")
print(kendall_tau)
cat("\nAll outputs saved.\n")
