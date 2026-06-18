#!/usr/bin/env Rscript
# ==============================================================================
# compute_all_metrics.R
# Comprehensive metrics computation for SVG rotation invariance benchmark
# ==============================================================================

# ==============================================================================
# SECTION 1: SETUP AND CONFIGURATION
# ==============================================================================

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(patchwork)

# Try to load ggVennDiagram; if unavailable, skip Venn diagrams
ggVennDiagram_available <- requireNamespace("ggVennDiagram", quietly = TRUE)
if (ggVennDiagram_available) {
  library(ggVennDiagram)
}

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

benchmark_root <- file.path(project_root, "src", "03_benchmark", "outputs")
metrics_root   <- file.path(project_root, "src", "04_metrics", "outputs")

tools <- list.files(benchmark_root)
tools <- setdiff(tools, "boostgp")

comparison_dir    <- file.path(metrics_root, "comparison")
comparison_figures <- file.path(comparison_dir, "figures")
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comparison_figures, recursive = TRUE, showWarnings = FALSE)

# --- Score configuration (for ranking/detection) ---
score_config <- list(
  sparkx    = list(col = "combinedPval", negate = TRUE),
  nnsvg     = list(col = "LR_stat",      negate = FALSE),
  spatialde = list(col = "FSV",          negate = FALSE),
  moransi   = list(col = "I",            negate = FALSE),
  smash     = list(col = "pval",         negate = TRUE)
)

# --- Adjusted p-value column for significance filtering ---
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
pairs <- list(c("0", "30"), c("0", "45"), c("0", "60"),
              c("30", "45"), c("30", "60"), c("45", "60"))
pair_labels <- sapply(pairs, function(p) paste0(p[1], "\u00b0 vs ", p[2], "\u00b0"))

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

# --- Helper to save both PNG and PDF ---
save_figure <- function(plot, filename_base, width, height, dpi = 300) {
  # Use cairo for proper Unicode support in PNG
  ggsave(paste0(filename_base, ".png"), plot, width = width, height = height,
         dpi = dpi, bg = "white", type = "cairo")
  ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height, bg = "white")
}


# ==============================================================================
# SECTION 2: HELPER FUNCTIONS
# ==============================================================================

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

compute_classification_metrics <- function(truth, predicted) {
  # Convert to numeric to prevent integer overflow in MCC denominator
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

  data.frame(
    TP = tp, TN = tn, FP = fp, FN = fn,
    Sensitivity = sensitivity, Specificity = specificity,
    Precision = precision, FDR = fdr, F1 = f1,
    MCC = mcc, BalancedAccuracy = balanced_acc,
    stringsAsFactors = FALSE
  )
}

# Classification tau: Spearman correlation between angle and metric value
# Measures whether classification performance degrades monotonically with rotation
compute_classification_tau <- function(angles_vec, metric_vec) {
  complete <- complete.cases(angles_vec, metric_vec)
  if (sum(complete) < 3) return(NA_real_)
  a <- angles_vec[complete]; m <- metric_vec[complete]
  if (length(unique(m)) < 2) return(NA_real_)
  cor(a, m, method = "spearman")
}

# Jaccard index for two sets
jaccard_index <- function(set_a, set_b) {
  inter <- length(intersect(set_a, set_b))
  uni <- length(union(set_a, set_b))
  if (uni == 0) return(NA_real_)
  inter / uni
}


# ==============================================================================
# SECTION 3: PER-METHOD ANALYSIS LOOP
# ==============================================================================

all_methods_metrics       <- list()
all_methods_rotation      <- list()
all_methods_rank_data     <- list()
all_methods_classification <- list()
all_methods_confusion      <- list()
all_methods_sig_sets       <- list()
all_methods_jaccard        <- list()
all_methods_sig_summary    <- list()

for (tool in tools) {
  cat("\n=====" , tool, "=====\n")

  benchmark_dir <- file.path(benchmark_root, tool)
  output_dir    <- file.path(metrics_root, tool)
  figures_dir   <- file.path(output_dir, "figures")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  cfg <- score_config[[tool]]
  pcol <- adj_pval_col[[tool]]
  if (is.null(cfg) || is.null(pcol)) {
    message("No config for ", tool, ", skipping")
    next
  }

  per_angle_scores <- list()
  per_angle_data   <- list()
  per_angle_class  <- list()
  per_angle_confusion <- list()
  sig_sets         <- list()

  for (angle in angles) {
    rds_file <- file.path(benchmark_dir, paste0("scdesign3_angle", angle, "_results.rds"))
    runtime_file <- file.path(benchmark_dir, paste0("scdesign3_angle", angle, "_runtime.csv"))
    if (!file.exists(rds_file) || !file.exists(runtime_file)) {
      message("Skipping angle ", angle, " -- missing outputs")
      next
    }

    results <- readRDS(rds_file)
    df <- results$res_mtest
    feature_names <- rownames(df)

    raw_score <- df[[cfg$col]]
    score <- if (cfg$negate) -raw_score else raw_score
    names(score) <- feature_names

    adj_p <- df[[pcol]]
    names(adj_p) <- feature_names

    parts <- strsplit(feature_names, "_")
    gene <- sapply(parts, `[`, 1)
    alpha <- as.numeric(sapply(parts, function(x) paste(x[-1], collapse = "_")))

    truth_any   <- as.integer(alpha > 0)
    truth_strong <- as.integer(alpha > 0.5)
    baseline_any <- sum(truth_any == 1L) / length(truth_any)
    baseline_strong <- sum(truth_strong == 1L) / length(truth_strong)

    auprc_any    <- compute_auprc(truth_any, score)
    auprc_strong <- compute_auprc(truth_strong, score)
    auroc_any    <- compute_auroc(truth_any, score)
    auroc_strong <- compute_auroc(truth_strong, score)

    tau_alpha_overall <- cor(alpha, score, method = "kendall", use = "complete.obs")
    tau_per_gene <- tapply(seq_along(gene), gene, function(idx) {
      if (length(idx) < 2) return(NA_real_)
      cor(alpha[idx], score[idx], method = "kendall", use = "complete.obs")
    })
    tau_alpha_gene <- mean(tau_per_gene, na.rm = TRUE)

    runtime <- read.csv(runtime_file)$elapsed_sec

    # --- Classification metrics at adj p < 0.05 ---
    predicted <- as.integer(adj_p < 0.05)
    predicted[is.na(predicted)] <- 0
    class_metrics <- compute_classification_metrics(truth_any, predicted)
    class_metrics$angle <- angle
    class_metrics$tool  <- tool
    per_angle_class[[as.character(angle)]] <- class_metrics

    confusion_df <- data.frame(
      tool = tool, angle = angle,
      TP = class_metrics$TP, TN = class_metrics$TN,
      FP = class_metrics$FP, FN = class_metrics$FN,
      stringsAsFactors = FALSE
    )
    per_angle_confusion[[as.character(angle)]] <- confusion_df

    # --- Significant gene sets for Venn/Jaccard ---
    sig_features <- feature_names[!is.na(adj_p) & adj_p < 0.05]
    sig_sets[[as.character(angle)]] <- sig_features

    per_angle_data[[as.character(angle)]] <- data.frame(
      gene = gene, alpha = alpha, feature = feature_names,
      score = score, rank = rank(-score, ties.method = "average"),
      adj_p = adj_p, significant = !is.na(adj_p) & adj_p < 0.05,
      angle = angle, stringsAsFactors = FALSE
    )

    per_angle_scores[[as.character(angle)]] <- data.frame(
      angle = angle,
      auprc_any = auprc_any, auprc_strong = auprc_strong,
      auprc_lift_any = auprc_any / baseline_any,
      auprc_lift_strong = auprc_strong / baseline_strong,
      auroc_any = auroc_any, auroc_strong = auroc_strong,
      tau_alpha_overall = tau_alpha_overall, tau_alpha_gene = tau_alpha_gene,
      runtime_sec = runtime,
      n_sig = length(sig_features),
      stringsAsFactors = FALSE
    )
  }

  # --- Combine per-angle results ---
  metrics <- do.call(rbind, per_angle_scores)
  metrics$tool <- tool
  all_methods_metrics[[tool]] <- metrics

  class_df <- do.call(rbind, per_angle_class)
  all_methods_classification[[tool]] <- class_df

  confusion_df <- do.call(rbind, per_angle_confusion)
  all_methods_confusion[[tool]] <- confusion_df

  all_methods_sig_sets[[tool]] <- sig_sets

  # --- Save per-method CSVs ---
  write.csv(metrics, file.path(output_dir, paste0(tool, "_metrics.csv")), row.names = FALSE)
  write.csv(class_df, file.path(output_dir, paste0(tool, "_classification.csv")), row.names = FALSE)
  write.csv(confusion_df, file.path(output_dir, paste0(tool, "_confusion_matrix.csv")), row.names = FALSE)

  angles_found <- names(per_angle_data)
  if (length(angles_found) < 2) next

  # --- Pairwise tau_rotation (all genes) ---
  pairwise_tau <- list()
  for (pair in pairs) {
    a <- pair[1]; b <- pair[2]
    if (!a %in% angles_found || !b %in% angles_found) next
    common <- intersect(per_angle_data[[a]]$feature, per_angle_data[[b]]$feature)
    da <- per_angle_data[[a]]; db <- per_angle_data[[b]]
    rownames(da) <- da$feature; rownames(db) <- db$feature
    da <- da[common, ]; db <- db[common, ]
    tau <- cor(da$score, db$score, method = "kendall", use = "complete.obs")
    pairwise_tau[[length(pairwise_tau) + 1]] <- data.frame(
      angle_a = as.numeric(a), angle_b = as.numeric(b),
      tau_rotation = tau, pair_label = paste0(a, "\u00b0 vs ", b, "\u00b0"),
      stringsAsFactors = FALSE
    )
  }

  # --- Pairwise tau_rotation on significant genes (at 0°) ---
  sig_at_0 <- sig_sets[["0"]]
  if (length(sig_at_0) > 2 && "0" %in% angles_found) {
    for (pair in pairs) {
      a <- pair[1]; b <- pair[2]
      if (a != "0" && b != "0") next  # Only pairs involving 0°
      if (!a %in% angles_found || !b %in% angles_found) next
      da <- per_angle_data[[a]]; db <- per_angle_data[[b]]
      rownames(da) <- da$feature; rownames(db) <- db$feature
      common <- intersect(da$feature, db$feature)
      sig_common <- intersect(sig_at_0, common)
      if (length(sig_common) < 3) next
      tau_sig <- cor(da[sig_common, "score"], db[sig_common, "score"],
                     method = "kendall", use = "complete.obs")
      # Append to pairwise_tau
      pairwise_tau[[length(pairwise_tau) + 1]] <- data.frame(
        angle_a = as.numeric(a), angle_b = as.numeric(b),
        tau_rotation = tau_sig,
        pair_label = paste0(a, "\u00b0 vs ", b, "\u00b0 (sig only)"),
        stringsAsFactors = FALSE
      )
    }
  }

  rotation_tau <- do.call(rbind, pairwise_tau)
  rotation_tau$tool <- tool
  all_methods_rotation[[tool]] <- rotation_tau
  write.csv(rotation_tau, file.path(output_dir, paste0(tool, "_rotation_tau.csv")), row.names = FALSE)

  # --- Set overlap metrics (Jaccard, consistency) ---
  if (length(sig_sets) >= 2) {
    jaccard_results <- list()
    for (pair in pairs) {
      a <- pair[1]; b <- pair[2]
      if (!a %in% names(sig_sets) || !b %in% names(sig_sets)) next
      jac <- jaccard_index(sig_sets[[a]], sig_sets[[b]])
      jaccard_results[[length(jaccard_results) + 1]] <- data.frame(
        angle_a = as.numeric(a), angle_b = as.numeric(b),
        jaccard = jac,
        pair_label = paste0(a, "\u00b0 vs ", b, "\u00b0"),
        stringsAsFactors = FALSE
      )
    }
    jaccard_df <- do.call(rbind, jaccard_results)
    jaccard_df$tool <- tool
    all_methods_jaccard[[tool]] <- jaccard_df
    write.csv(jaccard_df, file.path(output_dir, paste0(tool, "_jaccard.csv")), row.names = FALSE)

    # Consistency summary
    if (length(sig_sets) == 4) {
      consistent <- Reduce(intersect, sig_sets)
      union_all  <- Reduce(union, sig_sets)
      sig_summary <- data.frame(
        tool = tool,
        tool_label = toupper(tool),
        n_sig_0   = length(sig_sets[["0"]]),
        n_sig_30  = length(sig_sets[["30"]]),
        n_sig_45  = length(sig_sets[["45"]]),
        n_sig_60  = length(sig_sets[["60"]]),
        n_consistent = length(consistent),
        n_union      = length(union_all),
        consistency_rate = length(consistent) / max(1, length(union_all)),
        mean_jaccard = mean(jaccard_df$jaccard, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      all_methods_sig_summary[[tool]] <- sig_summary
      write.csv(sig_summary, file.path(output_dir, paste0(tool, "_sig_summary.csv")), row.names = FALSE)
    }
  }

  cat("Metrics:\n"); print(metrics)
  cat("\nRotation consistency:\n"); print(rotation_tau)
  cat("\nClassification metrics:\n"); print(class_df)

  all_data <- do.call(rbind, per_angle_data)
  all_data$angle <- factor(all_data$angle, levels = as.character(angles))

  if ("0" %in% angles_found && "60" %in% angles_found) {
    rank_0 <- per_angle_data[["0"]][, c("feature", "rank")]
    rank_60 <- per_angle_data[["60"]][, c("feature", "rank")]
    colnames(rank_0) <- c("feature", "rank_0")
    colnames(rank_60) <- c("feature", "rank_60")
    rank_merge <- merge(rank_0, rank_60, by = "feature")
    rank_merge$tool <- tool
    all_methods_rank_data[[tool]] <- rank_merge
  }

  # ==============================================================================
  # PER-METHOD PLOTS
  # ==============================================================================

  # --- AUC-ROC ---
  p1 <- ggplot(metrics, aes(x = angle)) +
    geom_line(aes(y = auroc_any, group = 1), color = method_colors[toupper(tool)], linewidth = 1) +
    geom_point(aes(y = auroc_any), color = method_colors[toupper(tool)], size = 3) +
    scale_y_continuous(limits = c(0.5, 1)) +
    labs(x = "Rotation angle (\u00b0)", y = "AUC-ROC",
         title = paste0(toupper(tool), ": AUC-ROC vs Rotation Angle"),
         subtitle = "Detection of spatially variable genes (\u03b1 > 0)") +
    theme_bioinformatics()
  save_figure(p1, file.path(figures_dir, paste0(tool, "_auroc")), width = 8, height = 6)

  # --- tau_alpha ---
  p2 <- ggplot(metrics, aes(x = angle, y = tau_alpha_gene)) +
    geom_line(group = 1, color = method_colors[toupper(tool)], linewidth = 1) +
    geom_point(color = method_colors[toupper(tool)], size = 3) +
    labs(x = "Rotation angle (\u00b0)", y = expression(tau[alpha] ~ "(per-gene mean)"),
         title = substitute(TOOL ~ ": " * tau[alpha] ~ "vs Rotation Angle", list(TOOL = toupper(tool))),
         subtitle = "Kendall \u03c4 between true \u03b1 and score, averaged per gene") +
    theme_bioinformatics()
  save_figure(p2, file.path(figures_dir, paste0(tool, "_tau")), width = 8, height = 6)

  # --- Pairwise scatter plots ---
  pairwise_df <- lapply(pairs, function(pair) {
    a <- pair[1]; b <- pair[2]
    if (!a %in% angles_found || !b %in% angles_found) return(NULL)
    da <- per_angle_data[[a]]; db <- per_angle_data[[b]]
    rownames(da) <- da$feature; rownames(db) <- db$feature
    common <- intersect(da$feature, db$feature)
    tau_val <- rotation_tau[rotation_tau$angle_a == as.numeric(a) &
                              rotation_tau$angle_b == as.numeric(b) &
                              !grepl("sig only", rotation_tau$pair_label), "tau_rotation"]
    if (length(tau_val) == 0) return(NULL)
    data.frame(feature = common, score_a = da[common, "score"], score_b = db[common, "score"],
               angle_a = a, angle_b = b, tau_rotation = tau_val,
               label = paste0(a, "\u00b0 vs ", b, "\u00b0\n\u03c4 = ", round(tau_val, 3)),
               stringsAsFactors = FALSE)
  })
  pairwise_df <- Filter(Negate(is.null), pairwise_df)
  if (length(pairwise_df) > 0) {
    all_pairs_df <- do.call(rbind, pairwise_df)
    all_pairs_df$label <- factor(all_pairs_df$label, levels = unique(all_pairs_df$label))
    lims_global <- range(c(all_pairs_df$score_a, all_pairs_df$score_b), na.rm = TRUE)

    for (pair in pairs) {
      a <- pair[1]; b <- pair[2]
      sub <- all_pairs_df[all_pairs_df$angle_a == a & all_pairs_df$angle_b == b, ]
      if (nrow(sub) == 0) next
      tau <- sub$tau_rotation[1]
      lims <- range(c(sub$score_a, sub$score_b), na.rm = TRUE)
      p <- ggplot(sub, aes(x = score_a, y = score_b)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
        geom_point(alpha = 0.4, size = 1.2, color = method_colors[toupper(tool)]) +
        coord_fixed(xlim = lims, ylim = lims) +
        labs(x = paste0("Statistic at ", a, "\u00b0"), y = paste0("Statistic at ", b, "\u00b0"),
             title = paste0(toupper(tool), ": Score Comparison ", a, "\u00b0 vs ", b, "\u00b0"),
             subtitle = paste0("Kendall \u03c4 = ", round(tau, 3))) +
        theme_bioinformatics()
      save_figure(p, file.path(figures_dir, paste0(tool, "_scatter_", a, "_vs_", b)),
                  width = 7, height = 7)
    }

    p_facet <- ggplot(all_pairs_df, aes(x = score_a, y = score_b)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(alpha = 0.3, size = 0.7, color = method_colors[toupper(tool)]) +
      facet_wrap(~ label, ncol = 3) +
      coord_fixed(xlim = lims_global, ylim = lims_global) +
      labs(x = "Statistic at first angle", y = "Statistic at second angle",
           title = paste0(toupper(tool), ": Pairwise Score Comparisons")) +
      theme_bioinformatics()
    save_figure(p_facet, file.path(figures_dir, paste0(tool, "_scatter_facet")),
                width = 12, height = 9)
  }

  # --- Classification metrics by angle ---
  if (nrow(class_df) > 0) {
    class_melt <- melt(class_df,
                       id.vars = c("tool", "angle"),
                       measure.vars = c("Sensitivity", "Specificity", "Precision",
                                        "FDR", "F1", "MCC", "BalancedAccuracy"),
                       variable.name = "Metric", value.name = "Value")
    class_melt$angle <- factor(class_melt$angle, levels = as.character(angles))

    p_class <- ggplot(class_melt, aes(x = angle, y = Value, group = Metric, color = Metric)) +
      geom_line(linewidth = 1) + geom_point(size = 2.5) +
      scale_color_brewer(palette = "Set1", name = "Metric") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Rotation angle (\u00b0)", y = "Metric value",
           title = paste0(toupper(tool), ": Classification Metrics vs Rotation Angle"),
           subtitle = "At adjusted p < 0.05 significance threshold") +
      theme_bioinformatics() +
      theme(legend.position = "right")
    save_figure(p_class, file.path(figures_dir, paste0(tool, "_classification_metrics")),
                width = 10, height = 6)

    # MCC specifically
    p_mcc <- ggplot(class_df, aes(x = factor(angle), y = MCC)) +
      geom_col(fill = method_colors[toupper(tool)], alpha = 0.8, color = "black", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.3f", MCC)), vjust = -0.5, size = 3.5) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Rotation angle (\u00b0)", y = "Matthews Correlation Coefficient (MCC)",
           title = paste0(toupper(tool), ": MCC vs Rotation Angle"),
           subtitle = "Higher MCC = better balanced classification performance") +
      theme_bioinformatics()
    save_figure(p_mcc, file.path(figures_dir, paste0(tool, "_mcc")), width = 8, height = 6)

    # Confusion matrix heatmap (at 0°)
    conf_0 <- confusion_df[confusion_df$angle == 0, ]
    if (nrow(conf_0) > 0) {
      conf_melt <- melt(conf_0, id.vars = c("tool", "angle"),
                        measure.vars = c("TP", "TN", "FP", "FN"),
                        variable.name = "Cell", value.name = "Count")
      p_conf <- ggplot(conf_melt, aes(x = Cell, y = tool, fill = Count)) +
        geom_tile(color = "white", linewidth = 1) +
        geom_text(aes(label = Count), size = 5, fontface = "bold") +
        scale_fill_gradient(low = "#FEE8C8", high = "#E6550D", name = "Count") +
        labs(x = "Confusion Matrix Cell", y = NULL,
             title = paste0(toupper(tool), ": Confusion Matrix at 0\u00b0"),
             subtitle = "Adjusted p < 0.05 significance threshold") +
        theme_bioinformatics() +
        theme(axis.text.y = element_text(face = "bold", size = 12))
      save_figure(p_conf, file.path(figures_dir, paste0(tool, "_confusion_matrix")),
                  width = 8, height = 4)
    }
  }

  # --- Venn diagram (4-way: 0, 30, 45, 60) ---
  if (ggVennDiagram_available && length(sig_sets) == 4) {
    venn_data <- list(
      `0`  = sig_sets[["0"]],
      `30` = sig_sets[["30"]],
      `45` = sig_sets[["45"]],
      `60` = sig_sets[["60"]]
    )
    p_venn <- ggVennDiagram(venn_data,
                            label = "count",
                            label_alpha = 0,
                            edge_size = 0.5) +
      scale_fill_gradient(low = "#F7F7F7", high = method_colors[toupper(tool)], name = "Count") +
      labs(title = paste0(toupper(tool), ": Significant SVG Overlap Across Rotations"),
           subtitle = "Adjusted p < 0.05 | Intersection = consistent across all angles") +
      theme_bioinformatics() +
      theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    save_figure(p_venn, file.path(figures_dir, paste0(tool, "_venn")), width = 10, height = 8)
  }

  # --- Number of significant genes per angle ---
  if (nrow(metrics) > 0) {
    p_nsig <- ggplot(metrics, aes(x = factor(angle), y = n_sig)) +
      geom_col(fill = method_colors[toupper(tool)], alpha = 0.8, color = "black", linewidth = 0.3) +
      geom_text(aes(label = n_sig), vjust = -0.5, size = 3.5) +
      labs(x = "Rotation angle (\u00b0)", y = "Number of significant genes",
           title = paste0(toupper(tool), ": Significant Genes per Angle"),
           subtitle = "Adjusted p < 0.05 threshold") +
      theme_bioinformatics()
    save_figure(p_nsig, file.path(figures_dir, paste0(tool, "_n_sig")), width = 8, height = 6)
  }

  cat("Plots saved to", figures_dir, "\n")
}


# ==============================================================================
# SECTION 4: CROSS-METHOD AGGREGATION
# ==============================================================================

cat("\n\n========== CROSS-METHOD COMPARISON ==========\n")

all_metrics <- do.call(rbind, all_methods_metrics)
all_metrics$tool_label <- toupper(all_metrics$tool)
write.csv(all_metrics, file.path(comparison_dir, "all_metrics.csv"), row.names = FALSE)

all_rotation <- do.call(rbind, all_methods_rotation)
all_rotation$tool_label <- toupper(all_rotation$tool)
all_rotation$pair_label <- factor(all_rotation$pair_label,
                                   levels = c(pair_labels, paste0(pair_labels, " (sig only)")))
write.csv(all_rotation, file.path(comparison_dir, "all_rotation_tau.csv"), row.names = FALSE)

all_classification <- do.call(rbind, all_methods_classification)
all_classification$tool_label <- toupper(all_classification$tool)
write.csv(all_classification, file.path(comparison_dir, "classification_metrics.csv"), row.names = FALSE)

all_confusion <- do.call(rbind, all_methods_confusion)
all_confusion$tool_label <- toupper(all_confusion$tool)
write.csv(all_confusion, file.path(comparison_dir, "confusion_matrix_summary.csv"), row.names = FALSE)

all_jaccard <- do.call(rbind, all_methods_jaccard)
all_jaccard$tool_label <- toupper(all_jaccard$tool)
write.csv(all_jaccard, file.path(comparison_dir, "all_jaccard.csv"), row.names = FALSE)

all_sig_summary <- do.call(rbind, all_methods_sig_summary)
write.csv(all_sig_summary, file.path(comparison_dir, "svg_consistency.csv"), row.names = FALSE)

# --- Rotation invariance summary (all genes) ---
rotation_main <- all_rotation[!grepl("sig only", all_rotation$pair_label), ]
mean_rotation <- aggregate(tau_rotation ~ tool_label, data = rotation_main,
                           FUN = function(x) c(mean = mean(x), min = min(x), max = max(x)))
mean_rotation <- do.call(data.frame, mean_rotation)
colnames(mean_rotation) <- c("tool_label", "mean_tau", "min_tau", "max_tau")
mean_rotation <- mean_rotation[order(mean_rotation$mean_tau, decreasing = TRUE), ]
write.csv(mean_rotation, file.path(comparison_dir, "rotation_invariance_summary.csv"), row.names = FALSE)

# --- Classification tau (Spearman correlation between angle and metric) ---
classification_tau_list <- list()
for (tool in unique(all_classification$tool)) {
  sub <- all_classification[all_classification$tool == tool, ]
  for (metric in c("Sensitivity", "Specificity", "Precision", "FDR", "F1", "MCC", "BalancedAccuracy")) {
    tau_val <- compute_classification_tau(sub$angle, sub[[metric]])
    classification_tau_list[[length(classification_tau_list) + 1]] <- data.frame(
      tool = tool, tool_label = toupper(tool),
      metric = metric, classification_tau = tau_val,
      stringsAsFactors = FALSE
    )
  }
}
classification_tau_df <- do.call(rbind, classification_tau_list)
write.csv(classification_tau_df, file.path(comparison_dir, "classification_tau.csv"), row.names = FALSE)


# ==============================================================================
# SECTION 5: CROSS-METHOD COMPARISON PLOTS
# ==============================================================================

# --- 5.1 Rotation invariance heatmap (all genes) ---
p_heat <- ggplot(rotation_main, aes(x = pair_label, y = tool_label, fill = tau_rotation)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", tau_rotation)), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#D73027", mid = "#FFFFCC", high = "#1A9850",
                       midpoint = 0.9, limits = c(0.5, 1), name = expression(tau[rotation])) +
  labs(x = "Angle pair", y = "Method",
       title = "Rotation Invariance: Score Consistency Across Angle Pairs",
       subtitle = "Higher \u03c4 = more rotation-invariant (scores stable after rotation)") +
  theme_bioinformatics() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_figure(p_heat, file.path(comparison_figures, "rotation_invariance_heatmap"), width = 10, height = 6)

# --- 5.2 Rotation invariance summary bar chart ---
p_rot_summary <- ggplot(mean_rotation, aes(x = reorder(tool_label, mean_tau), y = mean_tau)) +
  geom_col(fill = "#1A9850", alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = min_tau, ymax = max_tau), width = 0.2, color = "gray30", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.4f", mean_tau)), vjust = -0.5, size = 3.5, fontface = "bold") +
  labs(x = "Method", y = expression("Mean " * tau[rotation] ~ "(across all pairs)"),
       title = "Rotation Invariance Summary",
       subtitle = "Higher = more consistent scores across rotations (error bars: min/max pair)") +
  coord_flip(ylim = c(0.5, 1)) +
  theme_bioinformatics()
save_figure(p_rot_summary, file.path(comparison_figures, "rotation_invariance_summary"),
            width = 9, height = 6)

# --- 5.3 Detection vs Rotation Invariance trade-off ---
auroc_0 <- all_metrics[all_metrics$angle == 0, c("tool_label", "auroc_any")]
colnames(auroc_0)[2] <- "auroc_0"
tradeoff <- merge(auroc_0, mean_rotation, by = "tool_label")

p_tradeoff <- ggplot(tradeoff, aes(x = auroc_0, y = mean_tau, color = tool_label)) +
  geom_point(size = 5) +
  geom_text(aes(label = tool_label), vjust = -1.3, hjust = 0.5, size = 4, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = method_colors, guide = "none") +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "gray60", linewidth = 0.5) +
  labs(x = "AUC-ROC at 0\u00b0 (detection ability)",
       y = expression("Mean " * tau[rotation] ~ "(rotation invariance)"),
       title = "Detection vs Rotation Invariance Trade-off",
       subtitle = "Top-right = ideal: detects SVGs AND is rotation-invariant") +
  coord_cartesian(xlim = c(0.8, 1), ylim = c(0.7, 1)) +
  theme_bioinformatics()
save_figure(p_tradeoff, file.path(comparison_figures, "detection_vs_invariance"), width = 9, height = 8)

# --- 5.4 All methods AUC-ROC ---
p_all_auroc <- ggplot(all_metrics, aes(x = angle, y = auroc_any, color = tool_label, group = tool)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3) +
  scale_color_manual(values = method_colors, name = "Method") +
  labs(x = "Rotation angle (\u00b0)", y = "AUC-ROC",
       title = "AUC-ROC vs Rotation Angle \u2014 All Methods",
       subtitle = "Detection of spatially variable genes (\u03b1 > 0, baseline = 0.5)") +
  theme_bioinformatics() +
  theme(legend.position = "right")
save_figure(p_all_auroc, file.path(comparison_figures, "all_auroc"), width = 10, height = 6)

# --- 5.5 All methods tau_alpha ---
p_all_tau_gene <- ggplot(all_metrics, aes(x = angle, y = tau_alpha_gene, color = tool_label, group = tool)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3) +
  scale_color_manual(values = method_colors, name = "Method") +
  labs(x = "Rotation angle (\u00b0)", y = expression(tau[alpha] ~ "(per-gene mean)"),
       title = expression(tau[alpha] ~ "vs Rotation Angle \u2014 All Methods"),
       subtitle = "Concordance with true signal fraction, averaged per gene") +
  theme_bioinformatics() +
  theme(legend.position = "right")
save_figure(p_all_tau_gene, file.path(comparison_figures, "all_tau_per_gene"), width = 10, height = 6)

# --- 5.6 Runtime comparison ---
p_runtime <- ggplot(all_metrics, aes(x = factor(angle), y = runtime_sec, fill = tool_label)) +
  geom_col(position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = method_colors, name = "Method") +
  labs(x = "Rotation angle (\u00b0)", y = "Runtime (seconds)",
       title = "Runtime per Angle \u2014 All Methods") +
  theme_bioinformatics() +
  theme(legend.position = "right")
save_figure(p_runtime, file.path(comparison_figures, "all_runtime"), width = 10, height = 6)

# --- 5.7 Rank stability facet ---
if (length(all_methods_rank_data) > 0) {
  all_ranks <- do.call(rbind, all_methods_rank_data)
  all_ranks$tool_label <- toupper(all_ranks$tool)

  p_rank <- ggplot(all_ranks, aes(x = rank_0, y = rank_60)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(alpha = 0.3, size = 0.6, aes(color = tool_label)) +
    facet_wrap(~ tool_label, ncol = 3) +
    scale_color_manual(values = method_colors, guide = "none") +
    coord_fixed() +
    labs(x = "Rank at 0\u00b0", y = "Rank at 60\u00b0",
         title = "Rank Stability: 0\u00b0 vs 60\u00b0 \u2014 All Methods",
         subtitle = "Diagonal = perfect rotation invariance; scatter = rotation-variant") +
    theme_bioinformatics()
  save_figure(p_rank, file.path(comparison_figures, "rank_stability_facet"), width = 12, height = 9)
}


# ==============================================================================
# SECTION 6: NEW CLASSIFICATION COMPARISON PLOTS
# ==============================================================================

# --- 6.1 Classification metrics at 0° (grouped bar chart) ---
class_0 <- all_classification[all_classification$angle == 0, ]
class_0_melt <- melt(class_0,
                     id.vars = c("tool_label"),
                     measure.vars = c("Sensitivity", "Specificity", "Precision",
                                      "FDR", "F1", "MCC", "BalancedAccuracy"),
                     variable.name = "Metric", value.name = "Value")

p_class_0 <- ggplot(class_0_melt, aes(x = tool_label, y = Value, fill = tool_label)) +
  geom_col(color = "black", linewidth = 0.2) +
  facet_wrap(~ Metric, ncol = 4) +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Method", y = "Metric value",
       title = "Classification Metrics at 0\u00b0 \u2014 All Methods",
       subtitle = "Adjusted p < 0.05 significance threshold | Higher = better (except FDR)") +
  theme_bioinformatics() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_figure(p_class_0, file.path(comparison_figures, "classification_metrics_0deg"),
            width = 14, height = 8)

# --- 6.2 Classification metrics by angle (line plot, faceted by metric) ---
class_melt_all <- melt(all_classification,
                       id.vars = c("tool_label", "angle"),
                       measure.vars = c("Sensitivity", "Specificity", "Precision",
                                        "FDR", "F1", "MCC", "BalancedAccuracy"),
                       variable.name = "Metric", value.name = "Value")
class_melt_all$angle <- factor(class_melt_all$angle, levels = as.character(angles))

p_class_angle <- ggplot(class_melt_all, aes(x = angle, y = Value, color = tool_label, group = tool_label)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  facet_wrap(~ Metric, ncol = 4, scales = "free_y") +
  scale_color_manual(values = method_colors, name = "Method") +
  labs(x = "Rotation angle (\u00b0)", y = "Metric value",
       title = "Classification Metrics vs Rotation Angle \u2014 All Methods",
       subtitle = "Adjusted p < 0.05 | Stable lines = rotation-invariant classification") +
  theme_bioinformatics() +
  theme(legend.position = "bottom")
save_figure(p_class_angle, file.path(comparison_figures, "classification_metrics_by_angle"),
            width = 14, height = 10)

# --- 6.3 MCC comparison (single best summary metric) ---
p_mcc_comp <- ggplot(class_0, aes(x = reorder(tool_label, MCC), y = MCC, fill = tool_label)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", MCC)), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Method", y = "Matthews Correlation Coefficient (MCC)",
       title = "Classification Performance at 0\u00b0: MCC Comparison",
       subtitle = "MCC = +1 (perfect) to 0 (random) to -1 (anti-correlated) | Best metric for imbalanced data") +
  theme_bioinformatics()
save_figure(p_mcc_comp, file.path(comparison_figures, "mcc_comparison"), width = 9, height = 7)

# --- 6.4 Confusion matrix heatmap (all methods at 0°) ---
conf_0_all <- all_confusion[all_confusion$angle == 0, ]
conf_melt_all <- melt(conf_0_all, id.vars = c("tool_label", "angle"),
                      measure.vars = c("TP", "TN", "FP", "FN"),
                      variable.name = "Cell", value.name = "Count")

p_conf_all <- ggplot(conf_melt_all, aes(x = Cell, y = reorder(tool_label, Count), fill = Count)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Count), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "#FEE8C8", high = "#E6550D", name = "Count") +
  labs(x = "Confusion Matrix Cell", y = "Method",
       title = "Confusion Matrix at 0\u00b0 \u2014 All Methods",
       subtitle = "Adjusted p < 0.05 significance threshold | TN = 50 total (true negatives)") +
  theme_bioinformatics() +
  theme(axis.text.y = element_text(face = "bold", size = 11))
save_figure(p_conf_all, file.path(comparison_figures, "confusion_matrix_heatmap"),
            width = 10, height = 6)

# --- 6.5 Classification tau heatmap ---
class_tau_melt <- dcast(classification_tau_df, tool_label ~ metric, value.var = "classification_tau")
class_tau_melt2 <- melt(class_tau_melt, id.vars = "tool_label",
                        variable.name = "Metric", value.name = "Tau")

p_class_tau <- ggplot(class_tau_melt2, aes(x = Metric, y = tool_label, fill = Tau)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = ifelse(is.na(Tau), "NA", sprintf("%.3f", Tau))),
            size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#D73027", mid = "#FFFFCC", high = "#1A9850",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman \u03c1") +
  labs(x = "Classification Metric", y = "Method",
       title = "Classification Stability Across Rotations",
       subtitle = "Spearman \u03c1 between rotation angle and metric value | Negative = degrades with rotation") +
  theme_bioinformatics() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_figure(p_class_tau, file.path(comparison_figures, "classification_tau_heatmap"),
            width = 11, height = 6)


# ==============================================================================
# SECTION 7: NEW SET OVERLAP COMPARISON PLOTS
# ==============================================================================

# --- 7.1 Number of significant genes per angle ---
if (nrow(all_metrics) > 0) {
  p_nsig_all <- ggplot(all_metrics, aes(x = factor(angle), y = n_sig, fill = tool_label)) +
    geom_col(position = "dodge", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = method_colors, name = "Method") +
    labs(x = "Rotation angle (\u00b0)", y = "Number of significant genes",
         title = "Significant Genes per Angle \u2014 All Methods",
         subtitle = "Adjusted p < 0.05 threshold | Total genes = 1,050 (1,000 true positives, 50 true negatives)") +
    theme_bioinformatics() +
    theme(legend.position = "right")
  save_figure(p_nsig_all, file.path(comparison_figures, "n_sig_per_angle"), width = 10, height = 6)
}

# --- 7.2 SVG consistency summary ---
if (nrow(all_sig_summary) > 0) {
  p_consistency <- ggplot(all_sig_summary, aes(x = reorder(tool_label, consistency_rate), y = consistency_rate)) +
    geom_col(aes(fill = tool_label), color = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", 100 * consistency_rate)), vjust = -0.5,
              size = 4, fontface = "bold") +
    scale_fill_manual(values = method_colors, guide = "none") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(x = "Method", y = "Consistency rate (intersection / union)",
         title = "SVG Call Consistency Across All Rotations",
         subtitle = "Proportion of significant genes called at ALL angles (0\u00b0, 30\u00b0, 45\u00b0, 60\u00b0)") +
    theme_bioinformatics()
  save_figure(p_consistency, file.path(comparison_figures, "svg_consistency_summary"),
              width = 9, height = 7)

  # Also show n_consistent as absolute counts
  p_consistent_n <- ggplot(all_sig_summary, aes(x = reorder(tool_label, n_consistent), y = n_consistent)) +
    geom_col(aes(fill = tool_label), color = "black", linewidth = 0.3) +
    geom_text(aes(label = n_consistent), vjust = -0.5, size = 4, fontface = "bold") +
    scale_fill_manual(values = method_colors, guide = "none") +
    labs(x = "Method", y = "Consistent SVGs (intersection)",
         title = "Number of Consistent SVG Calls Across All Rotations",
         subtitle = "Genes called significant at ALL angles (0\u00b0, 30\u00b0, 45\u00b0, 60\u00b0)") +
    theme_bioinformatics()
  save_figure(p_consistent_n, file.path(comparison_figures, "svg_consistent_count"),
              width = 9, height = 7)
}

# --- 7.3 Jaccard heatmap ---
if (nrow(all_jaccard) > 0) {
  p_jaccard <- ggplot(all_jaccard, aes(x = pair_label, y = tool_label, fill = jaccard)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", jaccard)), size = 3.5, fontface = "bold") +
    scale_fill_gradient2(low = "#D73027", mid = "#FFFFCC", high = "#1A9850",
                         midpoint = 0.95, limits = c(0.9, 1), name = "Jaccard index") +
    labs(x = "Angle pair", y = "Method",
         title = "Set Overlap: Jaccard Index Across Angle Pairs",
         subtitle = "Jaccard = |intersection| / |union| of significant gene sets | Higher = more consistent calls") +
    theme_bioinformatics() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_figure(p_jaccard, file.path(comparison_figures, "jaccard_heatmap"), width = 10, height = 6)
}

# --- 7.4 Set vs Rank Consistency (THE key new plot) ---
if (nrow(all_sig_summary) > 0 && nrow(mean_rotation) > 0) {
  set_rank <- merge(all_sig_summary[, c("tool", "tool_label", "mean_jaccard", "consistency_rate")],
                    mean_rotation[, c("tool_label", "mean_tau")],
                    by = "tool_label")

  p_set_rank <- ggplot(set_rank, aes(x = mean_jaccard, y = mean_tau, color = tool_label)) +
    geom_point(size = 6) +
    geom_text(aes(label = tool_label), vjust = -1.3, hjust = 0.5, size = 4,
              fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = method_colors, guide = "none") +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "gray60", linewidth = 0.5) +
    geom_vline(xintercept = 0.95, linetype = "dashed", color = "gray60", linewidth = 0.5) +
    labs(x = "Mean Jaccard index (set-based consistency)",
         y = expression("Mean " * tau[rotation] ~ "(rank-based consistency)"),
         title = "Set Consistency vs Rank Consistency",
         subtitle = "Two dimensions of rotation invariance: are the SAME genes called? (x) Are they ranked the SAME way? (y)") +
    coord_cartesian(xlim = c(0.9, 1), ylim = c(0.5, 1)) +
    theme_bioinformatics()
  save_figure(p_set_rank, file.path(comparison_figures, "set_vs_rank_consistency"),
              width = 9, height = 8)
}


cat("\nAll comparison plots saved to", comparison_figures, "\n")
cat("Summary CSVs saved to", comparison_dir, "\n")
cat("\nRotation invariance ranking (rank-based):\n")
print(mean_rotation)

cat("\nSVG consistency ranking (set-based):\n")
if (nrow(all_sig_summary) > 0) {
  print(all_sig_summary[order(all_sig_summary$consistency_rate, decreasing = TRUE),
                        c("tool_label", "consistency_rate", "mean_jaccard")])
}

cat("\nClassification stability (Spearman \u03c1 between angle and MCC):\n")
print(classification_tau_df[classification_tau_df$metric == "MCC", c("tool_label", "classification_tau")])

cat("\nDone!\n")
