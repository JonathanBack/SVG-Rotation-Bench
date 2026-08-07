# ==============================================================================
# export_whole_dataset.R
# Exports the full unprocessed stxBrain slice (no scDesign3, no subsetting, no
# alpha blending) as counts.csv + location.csv to the "whole" mode output dir.
# Lets the benchmark methods run on the complete geneset so that Jaccard/Venn
# across rotations can be computed on the real data, in parallel to the
# simulated-data runs.
#
# Edit SLICES below to control which slices are processed.
# Output: src/01_simulation/outputs/scDesign3/{slice}/whole/data/{counts,location}.csv
# ==============================================================================

library(Seurat)
library(SeuratData)
library(SingleCellExperiment)

# --- Select which slices to run ---
SLICES <- c("anterior1", "anterior2", "posterior1", "posterior2")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

if (!"stxBrain" %in% InstalledData()) {
  message("Installing stxBrain SeuratData dataset...")
  InstallData("stxBrain")
}

for (slice in SLICES) {
  message("\n========== Exporting whole dataset: ", slice, " ==========\n")

  out_data_dir <- file.path(project_root, "src", "01_simulation", "outputs",
                            "scDesign3", slice, "whole", "data")
  dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)

  # Load stxBrain slice (full unprocessed gene set)
  seu <- LoadData("stxBrain", type = slice)

  # Spatial coordinates (lowres Visium grid) -> match simulation script convention
  coords <- GetTissueCoordinates(seu, scale = "lowres")
  df_loc <- data.frame(
    spatial1 = coords$x,
    spatial2 = coords$y,
    row.names = rownames(coords)
  )

  # Full counts matrix (genes x spots) - no subsetting, no filtering
  counts_mat <- GetAssayData(seu, assay = "Spatial", layer = "counts")

  message("  Loaded slice: ", slice, " (genes=", nrow(counts_mat),
          ", spots=", ncol(counts_mat), ")")

  # Export: feature names are plain gene symbols (no _alpha suffix)
  write.csv(as.data.frame(df_loc),
            file = file.path(out_data_dir, "location.csv"))
  write.csv(as.data.frame(as.matrix(counts_mat)),
            file = file.path(out_data_dir, "counts.csv"))

  message("  Done: ", slice, " -> ", out_data_dir)
}

message("\nAll whole-dataset exports complete.")