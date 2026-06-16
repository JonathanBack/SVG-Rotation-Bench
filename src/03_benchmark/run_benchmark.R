project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

source(file.path(project_root, "src", "03_benchmark", "benchmark_dispatch.R"))

spec <- benchmark_spec()
benchmark_output_dir <- file.path(project_root, spec$benchmark_dir)
results_manifest_file <- file.path(benchmark_output_dir, "results_manifest.csv")

dir.create(benchmark_output_dir, recursive = TRUE, showWarnings = FALSE)

results_manifest <- run_benchmark_jobs(spec = spec, project_root = project_root)
write.csv(results_manifest, results_manifest_file, row.names = FALSE)
