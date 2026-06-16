benchmark_spec <- function() {
  list(
    methods = c("sparkx", "smash", "spatialde", "nsvg"),
    angles_degrees = c(0, 30, 45, 60, 90),
    simulation_dir = file.path("src", "01_simulation", "outputs", "scDesign3"),
    rotation_dir = file.path("src", "02_rotation", "outputs"),
    benchmark_dir = file.path("src", "03_benchmark", "outputs")
  )
}

build_benchmark_jobs <- function(spec = benchmark_spec()) {
  expand.grid(
    method = spec$methods,
    angle_degrees = spec$angles_degrees,
    stringsAsFactors = FALSE
  )
}

normalize_result_path <- function(path, project_root = getwd()) {
  if (is.na(path) || length(path) == 0L || !nzchar(path)) {
    return(NA_character_)
  }

  if (grepl("^/", path)) {
    return(path)
  }

  file.path(project_root, path)
}

method_adapter_path <- function(method, project_root = getwd()) {
  file.path(project_root, "src", "03_benchmark", "methods", paste0(method, ".R"))
}

run_method_adapter <- function(method, counts_file, location_file, output_dir, angle_degrees, project_root = getwd()) {
  adapter_file <- method_adapter_path(method = method, project_root = project_root)

  if (!file.exists(adapter_file)) {
    stop(
      "No method adapter found for '", method, "'. Expected file: ", adapter_file,
      call. = FALSE
    )
  }

  source(adapter_file)

  runner_name <- paste0("run_", method)
  if (!exists(runner_name, mode = "function")) {
    stop(
      "Method adapter '", adapter_file, "' does not define function '", runner_name, "'.",
      call. = FALSE
    )
  }

  runner <- get(runner_name, mode = "function")
  runner(
    counts_file = counts_file,
    location_file = location_file,
    output_dir = output_dir,
    angle_degrees = angle_degrees,
    project_root = project_root
  )
}

run_benchmark_jobs <- function(spec = benchmark_spec(), project_root = getwd()) {
  jobs <- build_benchmark_jobs(spec)
  dir.create(file.path(project_root, spec$benchmark_dir), recursive = TRUE, showWarnings = FALSE)

  job_results <- vector("list", nrow(jobs))

  for (i in seq_len(nrow(jobs))) {
    method <- jobs$method[[i]]
    angle_degrees <- jobs$angle_degrees[[i]]
    location_file <- file.path(project_root, spec$rotation_dir, paste0("rotated_locations_", angle_degrees, ".csv"))
    counts_file <- file.path(project_root, spec$simulation_dir, "counts.csv")
    output_dir <- file.path(project_root, spec$benchmark_dir, method, paste0("angle_", angle_degrees))

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    job_results[[i]] <- tryCatch(
      {
        result <- run_method_adapter(
          method = method,
          counts_file = counts_file,
          location_file = location_file,
          output_dir = output_dir,
          angle_degrees = angle_degrees,
          project_root = project_root
        )

        data.frame(
          method = method,
          angle_degrees = angle_degrees,
          status = "ok",
          result_file = if (is.list(result) && !is.null(result$result_file)) {
            normalize_result_path(result$result_file, project_root = project_root)
          } else {
            NA_character_
          },
          runtime_sec = if (is.list(result) && !is.null(result$runtime_sec)) result$runtime_sec else NA_real_,
          peak_ram_gb = if (is.list(result) && !is.null(result$peak_ram_gb)) result$peak_ram_gb else NA_real_,
          stringsAsFactors = FALSE
        )
      },
      error = function(e) {
        data.frame(
          method = method,
          angle_degrees = angle_degrees,
          status = "error",
          error_message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  }

  do.call(rbind, job_results)
}
