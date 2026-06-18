library(anndata)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

boostgp_dir <- file.path(project_root, "src", "03_benchmark", "tools", "BOOST-GP")
anndata_dir <- file.path(project_root, "src", "02_rotation", "outputs", "anndata", "data")
output_dir <- file.path(project_root, "src", "03_benchmark", "outputs", "boostgp")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

owd <- setwd(boostgp_dir)
source("R/boost.gp.R")
setwd(owd)

angles_degrees <- c(0, 30, 45, 60)

for (angle in angles_degrees) {
  rds_file <- file.path(output_dir, paste0("scdesign3_angle", angle, "_results.rds"))
  runtime_file <- file.path(output_dir, paste0("scdesign3_angle", angle, "_runtime.csv"))

  if (file.exists(rds_file) && file.exists(runtime_file)) {
    message("Skipping angle ", angle, " -- outputs already exist")
    next
  }

  message("Running BOOST-GP for angle = ", angle)

  h5ad_file <- file.path(anndata_dir, paste0("scdesign3_angle", angle, ".h5ad"))
  adata <- read_h5ad(h5ad_file)

  counts <- as.matrix(adata$layers[["counts"]])
  colnames(counts) <- adata$var_names
  rownames(counts) <- adata$obs_names
  mode(counts) <- "integer"

  loc <- as.data.frame(adata$obsm[["spatial"]])
  rownames(loc) <- adata$obs_names
  colnames(loc) <- c("x", "y")

  t_start <- proc.time()
  result <- boost.gp(Y = counts, loc = loc, iter = 100, burn = 50)
  elapsed <- proc.time() - t_start

  result_list <- list(res_mtest = data.frame(
    adjustedPval = p.adjust(result$pval, method = "BH"),
    row.names = rownames(result)
  ))

  saveRDS(result_list, rds_file)
  write.csv(
    data.frame(angle = angle, elapsed_sec = unname(elapsed["elapsed"])),
    runtime_file,
    row.names = FALSE
  )

  message("Saved ", rds_file, " (", round(elapsed["elapsed"], 1), "s)")
}

message("BOOST-GP benchmark complete.")
