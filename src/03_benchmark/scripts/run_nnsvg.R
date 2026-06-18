library(anndata)
library(SpatialExperiment)
library(scran)
library(nnSVG)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

anndata_dir <- file.path(project_root, "src", "02_rotation", "outputs", "anndata", "data")
output_dir <- file.path(project_root, "src", "03_benchmark", "outputs", "nnsvg")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

angles_degrees <- c(0, 30, 45, 60)

for (angle in angles_degrees) {
  rds_file <- file.path(output_dir, paste0("scdesign3_angle", angle, "_results.rds"))
  runtime_file <- file.path(output_dir, paste0("scdesign3_angle", angle, "_runtime.csv"))

  if (file.exists(rds_file) && file.exists(runtime_file)) {
    message("Skipping angle ", angle, " -- outputs already exist")
    next
  }

  message("Running nnSVG for angle = ", angle)

  h5ad_file <- file.path(anndata_dir, paste0("scdesign3_angle", angle, ".h5ad"))
  adata <- read_h5ad(h5ad_file)

  counts <- t(as.matrix(adata$layers[["counts"]]))
  colnames(counts) <- adata$obs_names
  rownames(counts) <- adata$var_names

  loc <- as.data.frame(adata$obsm[["spatial"]])
  colnames(loc) <- c("x", "y")
  rownames(loc) <- colnames(counts)

  row_data <- adata$var
  row_data$gene_id <- rownames(row_data)
  row_data$feature_type <- "Gene Expression"

  spe <- SpatialExperiment(
    assays = list(counts = counts),
    rowData = row_data,
    colData = loc,
    spatialCoordsNames = c("x", "y")
  )

  spe <- computeLibraryFactors(spe)
  spe <- logNormCounts(spe)

  set.seed(2024)
  t_start <- proc.time()
  spe <- nnSVG(spe, n_threads = 5)
  elapsed <- proc.time() - t_start

  df <- rowData(spe)
  result <- list(res_mtest = data.frame(
    adjustedPval = df$padj,
    row.names = rownames(df)
  ))

  saveRDS(result, rds_file)
  write.csv(
    data.frame(angle = angle, elapsed_sec = unname(elapsed["elapsed"])),
    runtime_file,
    row.names = FALSE
  )

  message("Saved ", rds_file, " (", round(elapsed["elapsed"], 1), "s)")
}

message("nnSVG benchmark complete.")
