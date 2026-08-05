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

library(anndata)
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

    anndata_dir <- file.path(project_root, "src", "02_rotation", "outputs",
                             slice, mode, "anndata", "data")
    output_dir <- file.path(project_root, "src", "03_benchmark", "outputs",
                            slice, mode, "nnsvg")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    if (!dir.exists(anndata_dir)) {
      message("Skipping ", slice, "/", mode, " -- AnnData dir not found")
      next
    }
    message("\n=== nnSVG: ", slice, " / ", mode, " ===")

    for (angle in angles_degrees) {
      rds_file <- file.path(output_dir, paste0("scdesign3_angle", angle, "_results.rds"))
      runtime_file <- file.path(output_dir, paste0("scdesign3_angle", angle, "_runtime.csv"))

      # Skip if already computed
      if (file.exists(rds_file) && file.exists(runtime_file)) {
        message("  Skipping angle ", angle, " -- outputs already exist")
        next
      }

      h5ad_file <- file.path(anndata_dir, paste0("scdesign3_angle", angle, ".h5ad"))
      if (!file.exists(h5ad_file)) {
        message("  Missing AnnData for angle ", angle, ", skipping")
        next
      }
      message("  Running nnSVG for angle = ", angle)

      # --- Load pre-built AnnData ---
      adata <- read_h5ad(h5ad_file)

      # Transpose counts to genes x cells for SpatialExperiment
      counts <- t(as.matrix(adata$layers[["counts"]]))
      colnames(counts) <- adata$obs_names
      rownames(counts) <- adata$var_names

      # Extract spatial coordinates from .obsm
      loc <- as.data.frame(adata$obsm[["spatial"]])
      colnames(loc) <- c("x", "y")
      rownames(loc) <- colnames(counts)

      row_data <- adata$var
      row_data$gene_id <- rownames(row_data)
      row_data$feature_type <- "Gene Expression"

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