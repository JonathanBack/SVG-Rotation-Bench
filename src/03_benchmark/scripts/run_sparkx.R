# ==============================================================================
# run_sparkx.R
# Benchmarks the SPARK-X method for SVG detection across rotated datasets.
# For each (slice, mode, angle), loads the counts (shared across angles) and
# the rotated spatial locations, runs SPARK-X, and saves results.
#
# Edit SLICES and MODES below to control scope.
# Input:  src/01_simulation/outputs/scDesign3/{slice}/{mode}/data/counts.csv
#         src/02_rotation/outputs/{slice}/{mode}/locations/rotated_locations_{angle}.csv
# Output: src/03_benchmark/outputs/{slice}/{mode}/sparkx/scdesign3_angle{angle}_results.rds
#         src/03_benchmark/outputs/{slice}/{mode}/sparkx/scdesign3_angle{angle}_runtime.csv
# ==============================================================================

library(SPARK)

# --- Select which slices and modes to run ---
SLICES <- c("anterior1", "anterior2", "posterior1", "posterior2")
MODES  <- c("simulated", "whole")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

angles_degrees <- c(0, 30, 45, 60)

for (slice in SLICES) {
  for (mode in MODES) {

    counts_file <- file.path(project_root, "src", "01_simulation", "outputs",
                             "scDesign3", slice, mode, "data", "counts.csv")
    locations_dir <- file.path(project_root, "src", "02_rotation", "outputs",
                               slice, mode, "locations")
    output_dir <- file.path(project_root, "src", "03_benchmark", "outputs",
                            slice, mode, "sparkx")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    if (!file.exists(counts_file)) {
      message("Skipping ", slice, "/", mode, " -- counts.csv not found")
      next
    }
    if (!dir.exists(locations_dir)) {
      message("Skipping ", slice, "/", mode, " -- locations dir not found")
      next
    }
    message("\n=== SPARK-X: ", slice, " / ", mode, " ===")

    # Load count matrix (genes x cells); remove duplicate gene names if any
    counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
    counts <- as.matrix(counts)
    counts <- counts[!duplicated(rownames(counts)), ]

    # Filter zero-count genes (harmonize gene universe across all methods).
    # SPARK-X has an internal filter that drops zero-count genes; this makes
    # it explicit so all methods run on the identical gene set.
    counts <- counts[rowSums(counts) > 0, , drop = FALSE]

    for (angle in angles_degrees) {
      rds_file <- file.path(output_dir, paste0("scdesign3_angle", angle, "_results.rds"))
      runtime_file <- file.path(output_dir, paste0("scdesign3_angle", angle, "_runtime.csv"))

      if (file.exists(rds_file) && file.exists(runtime_file)) {
        message("  Skipping angle ", angle, " -- outputs already exist")
        next
      }

      location_file <- file.path(locations_dir, paste0("rotated_locations_", angle, ".csv"))
      if (!file.exists(location_file)) {
        message("  Missing locations for angle ", angle, ", skipping")
        next
      }
      message("  Running SPARK-X for angle = ", angle)

      locations <- read.csv(location_file, row.names = 1, check.names = FALSE)
      # Subset and order locations to match counts column order
      locs <- as.matrix(locations[colnames(counts), c("spatial1", "spatial2")])

      # --- Run SPARK-X with mixture model option ---
      t_start <- proc.time()
      result <- sparkx(counts, locs, numCores = 5, option = "mixture")
      elapsed <- proc.time() - t_start

      saveRDS(result, rds_file)
      write.csv(
        data.frame(angle = angle, elapsed_sec = unname(elapsed["elapsed"])),
        runtime_file,
        row.names = FALSE
      )

      message("    Saved ", rds_file, " (", round(elapsed["elapsed"], 1), "s)")
    }
  }
}

message("\nSPARK-X benchmark complete.")