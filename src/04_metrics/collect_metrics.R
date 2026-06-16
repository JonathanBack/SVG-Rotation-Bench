project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

source(file.path(project_root, "src", "04_metrics", "metrics_utils.R"))

benchmark_output_dir <- file.path(project_root, "src", "03_benchmark", "outputs")
metrics_output_dir <- file.path(project_root, "src", "04_metrics", "outputs")
results_manifest_file <- file.path(benchmark_output_dir, "results_manifest.csv")
metrics_file <- file.path(metrics_output_dir, "benchmark_metrics.csv")

dir.create(metrics_output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(results_manifest_file)) {
  stop(
    "Missing results manifest: ", results_manifest_file,
    call. = FALSE
  )
}

results_manifest <- read.csv(results_manifest_file, stringsAsFactors = FALSE)
metrics_table <- build_metric_summary(results_manifest)
write.csv(metrics_table, metrics_file, row.names = FALSE)
