# ==============================================================================
# run_nnsvg.R
# Benchmarks the nnSVG method for SVG detection across rotated datasets.
# For each (slice, mode, angle), loads the corresponding AnnData, builds a
# SpatialExperiment, normalizes, runs nnSVG, and saves results.
#
# Edit SLICES and MODES below to control scope.
# Input:  src/02_rotation/outputs/{slice}/{mode}/anndata/data/scdesign3_angle{angle}.h5ad
# Output: src/03_benchmark/outputs/{slice}/{mode}/nnsvg/scdesign3_angle{angle}_results.rds
#         src/03_benchmark/outputs/{slice}/{mode}/nnsvg/scdesign3_angle{angle}_runtime.csv
# ==============================================================================

library(SpatialExperiment)
library(scran)
library(nnSVG)

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
                            slice, mode, "nnsvg")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    if (!file.exists(counts_file)) {
      message("Skipping ", slice, "/", mode, " -- counts.csv not found")
      next
    }
    if (!dir.exists(locations_dir)) {
      message("Skipping ", slice, "/", mode, " -- locations dir not found")
      next
    }

    # Load counts once (genes x cells) — shared across all angles
    counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
    counts <- as.matrix(counts)

    # Minimal row metadata
    row_data <- data.frame(gene_name = rownames(counts))

    message("\n=== nnSVG: ", slice, " / ", mode, " ===")

    for (angle in angles_degrees) {
      rds_file <- file.path(output_dir, paste0("scdesign3_angle", angle, "_results.rds"))
      runtime_file <- file.path(output_dir, paste0("scdesign3_angle", angle, "_runtime.csv"))

      # Skip if already computed
      if (file.exists(rds_file) && file.exists(runtime_file)) {
        message("  Skipping angle ", angle, " -- outputs already exist")
        next
      }

      location_file <- file.path(locations_dir,
                                 paste0("rotated_locations_", angle, ".csv"))
      if (!file.exists(location_file)) {
        message("  Missing locations for angle ", angle, ", skipping")
        next
      }
      message("  Running nnSVG for angle = ", angle)

      # Load rotated spatial coordinates for this angle
      locations <- read.csv(location_file, row.names = 1, check.names = FALSE)

      # Coordinates for SpatialExperiment (renamed to x/y)
      loc <- as.data.frame(locations[, c("spatial1", "spatial2")])
      colnames(loc) <- c("x", "y")
      rownames(loc) <- colnames(counts)

      # Include cell_type in colData if present
      if ("cell_type" %in% names(locations)) {
        loc$cell_type <- locations$cell_type
      }

      # --- Build SpatialExperiment object ---
      spe <- SpatialExperiment(
        assays = list(counts = counts),
        rowData = row_data,
        colData = loc,
        spatialCoordsNames = c("x", "y")
      )

      # Normalize: compute size factors, then log-normalize
      spe <- computeLibraryFactors(spe)
      spe <- logNormCounts(spe)

      # Filter low-expressed and mitochondrial genes (nnSVG requirement)
      spe <- filter_genes(spe)

      # --- Run nnSVG ---
      set.seed(2024)
      t_start <- proc.time()
      spe <- nnSVG(spe, n_threads = 10)
      elapsed <- proc.time() - t_start

      # --- Wrap results for compatibility with downstream metric scripts ---
      df <- as.data.frame(rowData(spe))
      result <- list(res_mtest = df)

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

message("\nnnSVG benchmark complete.")