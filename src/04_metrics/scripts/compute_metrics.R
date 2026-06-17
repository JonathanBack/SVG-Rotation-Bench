library(ggplot2)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

benchmark_dir <- file.path(project_root, "src", "03_benchmark", "outputs", "sparkx")
output_dir <- file.path(project_root, "src", "04_metrics", "outputs")
figures_dir <- file.path(output_dir, "figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

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

angles <- c(0, 30, 45, 60)
per_angle_scores <- list()
per_angle_data <- list()

for (angle in angles) {
  rds_file <- file.path(benchmark_dir, paste0("scdesign3_angle", angle, "_results.rds"))
  runtime_file <- file.path(benchmark_dir, paste0("scdesign3_angle", angle, "_runtime.csv"))

  results <- readRDS(rds_file)
  pvals <- results$res_mtest$adjustedPval
  names(pvals) <- rownames(results$res_mtest)

  parts <- strsplit(names(pvals), "_")
  gene <- sapply(parts, `[`, 1)
  alpha <- as.numeric(sapply(parts, `[`, 2))

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
write.csv(metrics, file.path(output_dir, "sparkx_metrics.csv"), row.names = FALSE)

pairwise_tau <- list()
pairs <- list(
  c("0", "30"), c("0", "45"), c("0", "60"),
  c("30", "45"), c("30", "60"),
  c("45", "60")
)
for (pair in pairs) {
  a <- pair[1]
  b <- pair[2]
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
write.csv(rotation_tau, file.path(output_dir, "sparkx_rotation_tau.csv"), row.names = FALSE)

cat("Metrics:\n")
print(metrics)
cat("\nRotation consistency:\n")
print(rotation_tau)

all_data <- do.call(rbind, per_angle_data)
all_data$angle <- factor(all_data$angle, levels = as.character(angles))

p_metrics <- ggplot(metrics, aes(x = angle)) +
  geom_line(aes(y = auprc, group = 1), color = "#0072B2", linewidth = 1) +
  geom_point(aes(y = auprc), color = "#0072B2", size = 3) +
  geom_line(aes(y = tau_alpha + 1, group = 1), color = "#D55E00", linewidth = 1) +
  geom_point(aes(y = tau_alpha + 1), color = "#D55E00", size = 3) +
  scale_y_continuous(
    name = "auPRC",
    limits = c(0, 2),
    breaks = seq(0, 1, 0.25),
    sec.axis = sec_axis(transform = ~ . - 1, name = expression(tau[alpha]),
                         breaks = seq(-1, 1, 0.5))
  ) +
  labs(x = "Rotation angle", title = "SPARK-X: auPRC and \u03c4\u03b1 vs. Rotation Angle") +
  theme_minimal()

ggsave(file.path(figures_dir, "sparkx_metrics.png"), p_metrics, width = 6, height = 4, dpi = 150)

p_density <- ggplot(all_data, aes(x = score, fill = angle, color = angle)) +
  geom_density(alpha = 0.25, linewidth = 0.5) +
  labs(
    x = "score (-adjusted P-value)",
    y = "Density",
    title = "SPARK-X: Score Distribution by Rotation Angle"
  ) +
  theme_minimal()

ggsave(file.path(figures_dir, "sparkx_score_density.png"), p_density, width = 6, height = 4, dpi = 150)

pair_0_60 <- all_data[all_data$angle %in% c("0", "60"), ]
pair_0 <- pair_0_60[pair_0_60$angle == "0", ]
pair_60 <- pair_0_60[pair_0_60$angle == "60", ]
rownames(pair_0) <- pair_0$feature
rownames(pair_60) <- pair_60$feature
common_features <- intersect(pair_0$feature, pair_60$feature)

df_scatter <- data.frame(
  feature = common_features,
  score_0 = pair_0[common_features, "score"],
  score_60 = pair_60[common_features, "score"],
  stringsAsFactors = FALSE
)

lims <- range(c(df_scatter$score_0, df_scatter$score_60))

p_scatter <- ggplot(df_scatter, aes(x = score_0, y = score_60)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.5, size = 1.5, color = "#0072B2") +
  coord_fixed(xlim = lims, ylim = lims) +
  labs(
    x = "score at 0\u00b0",
    y = "score at 60\u00b0",
    title = "SPARK-X: Score Comparison 0\u00b0 vs 60\u00b0"
  ) +
  theme_minimal()

ggsave(file.path(figures_dir, "sparkx_pairwise_scatter.png"), p_scatter, width = 5, height = 5, dpi = 150)

cat("\nPlots saved to", figures_dir, "\n")
