library(ggplot2)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

benchmark_root <- file.path(project_root, "src", "03_benchmark", "outputs")
metrics_root <- file.path(project_root, "src", "04_metrics", "outputs")

tools <- list.files(benchmark_root)
tools <- setdiff(tools, "boostgp")

comparison_dir <- file.path(metrics_root, "comparison")
comparison_figures <- file.path(comparison_dir, "figures")
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comparison_figures, recursive = TRUE, showWarnings = FALSE)

compute_auprc <- function(truth, score) {
  truth <- as.integer(truth)
  score <- as.numeric(score)

  if (length(truth) == 0L || all(is.na(truth)) || length(unique(truth)) < 2L) {
    return(NA_real_)
  }

  valid_idx <- which(!is.na(truth) & !is.na(score))
  truth <- truth[valid_idx]
  score <- score[valid_idx]

  if (!any(truth == 1L) || !any(truth == 0L)) {
    return(NA_real_)
  }

  ordering <- order(score, decreasing = TRUE)
  truth <- truth[ordering]

  cumulative_positives <- cumsum(truth == 1L)
  precision <- cumulative_positives / seq_along(truth)
  recall <- cumulative_positives / sum(truth == 1L)
  recall_previous <- c(0, recall[-length(recall)])

  sum((recall - recall_previous) * precision)
}

all_methods_metrics <- list()
all_methods_rotation <- list()

angles <- c(0, 30, 45, 60)

for (tool in tools) {
  cat("\n=====", tool, "=====\n")

  benchmark_dir <- file.path(benchmark_root, tool)
  output_dir <- file.path(metrics_root, tool)
  figures_dir <- file.path(output_dir, "figures")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  per_angle_scores <- list()
  per_angle_data <- list()

  for (angle in angles) {
    rds_file <- file.path(benchmark_dir, paste0("scdesign3_angle", angle, "_results.rds"))
    runtime_file <- file.path(benchmark_dir, paste0("scdesign3_angle", angle, "_runtime.csv"))

    if (!file.exists(rds_file) || !file.exists(runtime_file)) {
      message("Skipping angle ", angle, " -- missing outputs")
      next
    }

    results <- readRDS(rds_file)
    pvals <- results$res_mtest$adjustedPval
    names(pvals) <- rownames(results$res_mtest)

    parts <- strsplit(names(pvals), "_")
    gene <- sapply(parts, `[`, 1)
    alpha <- as.numeric(sapply(parts, function(x) paste(x[-1], collapse = "_")))

    score <- -pvals
    truth <- as.integer(alpha > 0)

    auprc <- compute_auprc(truth, score)
    tau_alpha <- cor(alpha, score, method = "kendall", use = "complete.obs")
    runtime <- read.csv(runtime_file)$elapsed_sec

    per_angle_data[[as.character(angle)]] <- data.frame(
      gene = gene,
      alpha = alpha,
      feature = names(pvals),
      score = score,
      angle = angle,
      stringsAsFactors = FALSE
    )

    per_angle_scores[[as.character(angle)]] <- data.frame(
      angle = angle,
      auprc = auprc,
      tau_alpha = tau_alpha,
      runtime_sec = runtime,
      stringsAsFactors = FALSE
    )
  }

  metrics <- do.call(rbind, per_angle_scores)
  metrics$tool <- tool

  all_methods_metrics[[tool]] <- metrics

  write.csv(metrics, file.path(output_dir, paste0(tool, "_metrics.csv")), row.names = FALSE)

  angles_found <- names(per_angle_data)
  if (length(angles_found) < 2) next

  pairs <- list(
    c("0", "30"), c("0", "45"), c("0", "60"),
    c("30", "45"), c("30", "60"),
    c("45", "60")
  )

  pairwise_tau <- list()
  for (pair in pairs) {
    a <- pair[1]; b <- pair[2]
    if (!a %in% angles_found || !b %in% angles_found) next
    tau <- cor(
      per_angle_data[[a]]$score,
      per_angle_data[[b]]$score,
      method = "kendall", use = "complete.obs"
    )
    pairwise_tau[[length(pairwise_tau) + 1]] <- data.frame(
      angle_a = as.numeric(a),
      angle_b = as.numeric(b),
      tau_rotation = tau,
      stringsAsFactors = FALSE
    )
  }
  rotation_tau <- do.call(rbind, pairwise_tau)
  rotation_tau$tool <- tool

  all_methods_rotation[[tool]] <- rotation_tau

  write.csv(rotation_tau, file.path(output_dir, paste0(tool, "_rotation_tau.csv")), row.names = FALSE)

  cat("Metrics:\n")
  print(metrics)
  cat("\nRotation consistency:\n")
  print(rotation_tau)

  all_data <- do.call(rbind, per_angle_data)
  all_data$angle <- factor(all_data$angle, levels = as.character(angles))

  p_auprc <- ggplot(metrics, aes(x = angle, y = auprc)) +
    geom_line(group = 1, color = "#0072B2", linewidth = 1) +
    geom_point(color = "#0072B2", size = 3) +
    labs(
      x = "Rotation angle (degrees)",
      y = "auPRC",
      title = paste0(toupper(tool), ": auPRC vs Rotation Angle"),
      subtitle = "Area under Precision-Recall curve for detecting rotated SV genes"
    ) +
    theme_minimal()

  ggsave(file.path(figures_dir, paste0(tool, "_auprc.png")),
         p_auprc, width = 6, height = 4, dpi = 150)

  p_tau <- ggplot(metrics, aes(x = angle, y = tau_alpha)) +
    geom_line(group = 1, color = "#D55E00", linewidth = 1) +
    geom_point(color = "#D55E00", size = 3) +
    labs(
      x = "Rotation angle (degrees)",
      y = expression(tau[alpha]),
      title = substitute(TOOL ~ ": " * tau[alpha] ~ "vs Rotation Angle",
                         list(TOOL = toupper(tool))),
      subtitle = expression(Kendall ~ tau ~ "between true signal fraction" ~ alpha ~ "and score")
    ) +
    theme_minimal()

  ggsave(file.path(figures_dir, paste0(tool, "_tau.png")),
         p_tau, width = 6, height = 4, dpi = 150)

  pairwise_df <- lapply(pairs, function(pair) {
    a <- pair[1]; b <- pair[2]
    if (!a %in% angles_found || !b %in% angles_found) return(NULL)
    da <- per_angle_data[[a]]
    db <- per_angle_data[[b]]
    rownames(da) <- da$feature
    rownames(db) <- db$feature
    common <- intersect(da$feature, db$feature)
    tau <- rotation_tau[rotation_tau$angle_a == as.numeric(a) &
                        rotation_tau$angle_b == as.numeric(b), "tau_rotation"]
    if (length(tau) == 0) return(NULL)
    data.frame(
      feature = common,
      score_a = da[common, "score"],
      score_b = db[common, "score"],
      angle_a = a,
      angle_b = b,
      tau_rotation = tau,
      label = paste0(a, "\u00b0 vs ", b, "\u00b0\n\u03c4 = ", round(tau, 3)),
      stringsAsFactors = FALSE
    )
  })
  pairwise_df <- Filter(Negate(is.null), pairwise_df)
  if (length(pairwise_df) == 0) next

  all_pairs_df <- do.call(rbind, pairwise_df)
  all_pairs_df$label <- factor(all_pairs_df$label, levels = unique(all_pairs_df$label))

  lims_global <- range(c(all_pairs_df$score_a, all_pairs_df$score_b))

  for (pair in pairs) {
    a <- pair[1]; b <- pair[2]
    sub <- all_pairs_df[all_pairs_df$angle_a == a & all_pairs_df$angle_b == b, ]
    if (nrow(sub) == 0) next
    tau <- sub$tau_rotation[1]
    lims <- range(c(sub$score_a, sub$score_b))

    p <- ggplot(sub, aes(x = score_a, y = score_b)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(alpha = 0.5, size = 1.5, color = "#0072B2") +
      coord_fixed(xlim = lims, ylim = lims) +
      labs(
        x = paste0("Score (\u2212adjusted p-value) at ", a, "\u00b0"),
        y = paste0("Score (\u2212adjusted p-value) at ", b, "\u00b0"),
        title = paste0(toupper(tool), ": Score Comparison ", a, "\u00b0 vs ", b, "\u00b0"),
        subtitle = paste0("Kendall \u03c4 = ", round(tau, 3))
      ) +
      theme_minimal()

    ggsave(file.path(figures_dir, paste0(tool, "_scatter_", a, "_vs_", b, ".png")),
           p, width = 5, height = 5, dpi = 150)
  }

  p_facet <- ggplot(all_pairs_df, aes(x = score_a, y = score_b)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(alpha = 0.5, size = 0.8, color = "#0072B2") +
    facet_wrap(~ label, ncol = 3) +
    coord_fixed(xlim = lims_global, ylim = lims_global) +
    labs(
      x = "Score (\u2212adjusted p-value) at first angle",
      y = "Score (\u2212adjusted p-value) at second angle",
      title = paste0(toupper(tool), ": Pairwise Score Comparisons Across Rotation Angles")
    ) +
    theme_minimal()

  ggsave(file.path(figures_dir, paste0(tool, "_scatter_facet.png")),
         p_facet, width = 10, height = 7, dpi = 150)

  cat("Plots saved to", figures_dir, "\n")
}

cat("\n\n========== CROSS-METHOD COMPARISON ==========\n")

all_metrics <- do.call(rbind, all_methods_metrics)
all_metrics$tool_label <- toupper(all_metrics$tool)
write.csv(all_metrics, file.path(comparison_dir, "all_metrics.csv"), row.names = FALSE)

all_rotation <- do.call(rbind, all_methods_rotation)
all_rotation$tool_label <- toupper(all_rotation$tool)
write.csv(all_rotation, file.path(comparison_dir, "all_rotation_tau.csv"), row.names = FALSE)

method_colors <- c(
  SPARKX = "#0072B2", NNSVG = "#D55E00", SPATIALDE = "#009E73",
  MORANSI = "#CC79A7", SMASH = "#F0E442"
)

p_all_auprc <- ggplot(all_metrics, aes(x = angle, y = auprc, color = tool_label, group = tool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = method_colors, name = "Method") +
  labs(
    x = "Rotation angle (degrees)",
    y = "auPRC",
    title = "auPRC vs Rotation Angle — All Methods",
    subtitle = "Higher = better detection of spatially variable genes"
  ) +
  theme_minimal()

ggsave(file.path(comparison_figures, "all_auprc.png"),
       p_all_auprc, width = 7, height = 4.5, dpi = 150)

p_all_tau <- ggplot(all_metrics, aes(x = angle, y = tau_alpha, color = tool_label, group = tool)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = method_colors, name = "Method") +
  labs(
    x = "Rotation angle (degrees)",
    y = expression(tau[alpha]),
    title = expression(tau[alpha] ~ "vs Rotation Angle — All Methods"),
    subtitle = "Higher = better concordance with true signal fraction"
  ) +
  theme_minimal()

ggsave(file.path(comparison_figures, "all_tau.png"),
       p_all_tau, width = 7, height = 4.5, dpi = 150)

p_runtime <- ggplot(all_metrics, aes(x = factor(angle), y = runtime_sec, fill = tool_label)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = method_colors, name = "Method") +
  labs(
    x = "Rotation angle (degrees)",
    y = "Runtime (seconds)",
    title = "Runtime per Angle — All Methods"
  ) +
  theme_minimal()

ggsave(file.path(comparison_figures, "all_runtime.png"),
       p_runtime, width = 7, height = 4.5, dpi = 150)

rot_0_60 <- all_rotation[all_rotation$angle_a == 0 & all_rotation$angle_b == 60, ]
p_rot_tau <- ggplot(rot_0_60, aes(x = tool_label, y = tau_rotation, fill = tool_label)) +
  geom_col() +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(
    x = "Method",
    y = expression(tau[rotation]),
    title = "Rotation Consistency: 0\u00b0 vs 60\u00b0",
    subtitle = "Higher = scores more consistent after rotation"
  ) +
  theme_minimal()

ggsave(file.path(comparison_figures, "all_rotation_tau.png"),
       p_rot_tau, width = 6, height = 4.5, dpi = 150)

cat("\nAll comparison plots saved to", comparison_figures, "\n")
cat("Summary CSVs saved to", comparison_dir, "\n")
