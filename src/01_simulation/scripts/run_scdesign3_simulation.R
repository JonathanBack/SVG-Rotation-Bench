# ==============================================================================
# run_scdesign3_simulation.R
# Fits scDesign3 to each stxBrain slice (anterior1/2, posterior1/2), then
# generates the alpha-sweep simulated counts (alpha = 0, 0.05, ..., 1.0) by
# blending the fitted spatial mean with a shuffled (non-spatial) mean matrix.
#
# Per slice:
#   1. Load stxBrain slice via SeuratData; SCTransform + Moran's I (top 200)
#   2. First pass: construct_data -> fit_marginal (GP k=500, NB) -> fit_copula
#   3. Select top 50 genes by deviance explained
#   4. Second pass on the 50 genes
#   5. Sanity plot: Mbp real vs simulated at 100% signal
#   6. Alpha sweep -> stacked counts matrix with "gene_alpha" rownames
#   7. Save counts.csv + location.csv
#
# Edit SLICES below to control which slices are processed.
# Output: src/01_simulation/outputs/scDesign3/{slice}/simulated/data/{counts,location}.csv
# ==============================================================================

library(Seurat)
library(SeuratData)
library(SingleCellExperiment)
library(scDesign3)
library(scales)
library(ggplot2)
library(cowplot)
library(dplyr)

# --- Select which slices to run (edit these vectors to control scope) ---
SLICES <- c("anterior1", "anterior2", "posterior1", "posterior2")

# --- Setup: resolve project root ---
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Ensure stxBrain dataset is installed
if (!"stxBrain" %in% InstalledData()) {
  message("Installing stxBrain SeuratData dataset...")
  InstallData("stxBrain")
}

# ==============================================================================
# HELPERS
# ==============================================================================

# Save a plot as both PNG and PDF
save_plot <- function(plot_object, file_stub, width = 6, height = 5) {
  ggsave(
    filename = glue::glue("{file_stub}.png"),
    plot = plot_object, width = width, height = height, dpi = 300
  )
  ggsave(
    filename = glue::glue("{file_stub}.pdf"),
    plot = plot_object, width = width, height = height
  )
}

# Plot spatial expression of a single gene (rescaled log1p expression)
plot_exp <- function(sce, gene, pt_size = 1) {
  df_loc <- colData(sce)[, c("spatial1", "spatial2")]
  df_exp <- as.data.frame(counts(sce)[gene, ])
  colnames(df_exp) <- c("exp")
  df_exp$exp <- rescale(log1p(df_exp$exp))
  df <- cbind(df_exp, df_loc)

  p <- ggplot(data = df, aes(x = .data$spatial1, y = .data$spatial2)) +
    geom_point(aes(x = .data$spatial1, y = .data$spatial2, color = .data$exp),
               size = pt_size) +
    scale_colour_gradientn(colors = viridis_pal(option = "magma")(10),
                           limits = c(0, 1)) +
    theme_cowplot() +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    ggtitle(gene)
  return(p)
}

# ==============================================================================
# MAIN LOOP OVER SLICES
# ==============================================================================

for (slice in SLICES) {
  message("\n========== Processing slice: ", slice, " ==========\n")

  sim_data_dir <- file.path(project_root, "src", "01_simulation", "outputs",
                            "scDesign3", slice, "simulated", "data")
  sim_fig_dir  <- file.path(project_root, "src", "01_simulation", "outputs",
                            "scDesign3", slice, "simulated", "figures")
  dir.create(sim_data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(sim_fig_dir,  recursive = TRUE, showWarnings = FALSE)

  # ----------------------------------------------------------------------------
  # STEP 1: Load stxBrain, run Seurat-canonical preprocessing (SCTransform +
  # Moran's I) and convert the subsetted object to SingleCellExperiment
  # ----------------------------------------------------------------------------

  seu <- LoadData("stxBrain", type = slice)

  # Canonical Visium preprocessing: SCTransform performs variance-stabilizing
  # normalization, creates the "SCT" assay and sets VariableFeatures.
  seu <- SCTransform(seu, assay = "Spatial", verbose = FALSE)

  # Moran's I on the top 1000 variable features. Coordinates are read from the
  # Visium spatial image attached to the Seurat object automatically.
  seu <- FindSpatiallyVariableFeatures(
    seu,
    features = VariableFeatures(seu)[1:1000],
    selection.method = "moransi"
  )
  top.features <- head(SpatiallyVariableFeatures(seu, method = "moransi"), 200)

  # Subset the Seurat object to the top 200 SVGs
  seu <- seu[top.features, ]

  # Convert the subsetted Seurat object to SCE. Extract RAW counts from the
  # original "Spatial" assay (not the "SCT" assay) since scDesign3 needs raw
  # counts for construct_data(assay_use = "counts").
  coords <- GetTissueCoordinates(seu, scale = "lowres")
  # imagecol -> spatial1, imagerow -> spatial2 (match Visium orientation)
  df_loc <- data.frame(
    spatial1 = coords$imagecol,
    spatial2 = coords$imagerow,
    row.names = rownames(coords)
  )
  counts_mat <- GetAssayData(seu, assay = "Spatial", layer = "counts")

  ref_sce <- SingleCellExperiment(
    list(counts = as.matrix(counts_mat)),
    colData = df_loc
  )

  message("  Loaded slice: ", slice, " (top-200 SVGs, spots=", ncol(ref_sce), ")")

  # ----------------------------------------------------------------------------
  # STEP 2: First pass scDesign3 fit (top 200 genes)
  # ----------------------------------------------------------------------------

  set.seed(2024)

  ref_data <- construct_data(
    sce = ref_sce,
    assay_use = "counts",
    celltype = NULL,
    pseudotime = NULL,
    spatial = c("spatial1", "spatial2"),
    other_covariates = NULL,
    corr_by = "1"
  )

  ref_marginal <- fit_marginal(
    data = ref_data,
    predictor = "gene",
    mu_formula = "s(spatial1, spatial2, bs = 'gp', k = 500)",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 2,
    usebam = FALSE,
    trace = TRUE
  )

  ref_copula <- fit_copula(
    sce = ref_sce,
    assay_use = "counts",
    marginal_list = ref_marginal,
    family_use = "nb",
    copula = "gaussian",
    n_cores = 2,
    input_data = ref_data$dat
  )

  ref_para <- extract_para(
    sce = ref_sce,
    marginal_list = ref_marginal,
    n_cores = 5,
    family_use = "nb",
    new_covariate = ref_data$newCovariate,
    data = ref_data$dat
  )

  # ----------------------------------------------------------------------------
  # STEP 3: Select top 50 genes by deviance explained and subset
  # ----------------------------------------------------------------------------

  dev_explain <- sapply(ref_marginal, function(x) {
    sum <- summary(x$fit)
    return(sum$dev.expl)
  })
  dev_ordered <- order(dev_explain, decreasing = TRUE)
  num_de <- 50
  ordered <- dev_explain[dev_ordered]
  sel_genes <- names(ordered)[1:num_de]

  ref_sce <- ref_sce[sel_genes, ]
  message("  Subselected to top 50 genes by dev.expl")

  # ----------------------------------------------------------------------------
  # STEP 4: Second pass scDesign3 fit on top 50 genes
  # ----------------------------------------------------------------------------

  ref_data <- construct_data(
    sce = ref_sce,
    assay_use = "counts",
    celltype = NULL,
    pseudotime = NULL,
    spatial = c("spatial1", "spatial2"),
    other_covariates = NULL,
    corr_by = "1"
  )

  ref_marginal <- fit_marginal(
    data = ref_data,
    predictor = "gene",
    mu_formula = "s(spatial1, spatial2, bs = 'gp', k = 500)",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 2,
    usebam = FALSE,
    trace = TRUE
  )

  ref_copula <- fit_copula(
    sce = ref_sce,
    assay_use = "counts",
    marginal_list = ref_marginal,
    family_use = "nb",
    copula = "gaussian",
    n_cores = 2,
    input_data = ref_data$dat
  )

  ref_para <- extract_para(
    sce = ref_sce,
    marginal_list = ref_marginal,
    n_cores = 5,
    family_use = "nb",
    new_covariate = ref_data$newCovariate,
    data = ref_data$dat
  )

  # ----------------------------------------------------------------------------
  # STEP 5: Sanity plot - Mbp real vs simulated at 100% signal
  # ----------------------------------------------------------------------------

  sim_count_full <- simu_new(
    sce = ref_sce,
    mean_mat = ref_para$mean_mat,
    sigma_mat = ref_para$sigma_mat,
    zero_mat = ref_para$zero_mat,
    quantile_mat = NULL,
    copula_list = ref_copula$copula_list,
    n_cores = 5,
    family_use = "nb",
    input_data = ref_data$dat,
    new_covariate = ref_data$newCovariate,
    important_feature = rep(TRUE, dim(ref_sce)[1]),
    filtered_gene = NULL
  )

  sim_sce <- SingleCellExperiment(list(counts = sim_count_full),
                                  colData = ref_data$newCovariate)

  sanity_gene <- if ("Mbp" %in% rownames(ref_sce)) "Mbp" else rownames(ref_sce)[1]
  message("  Generating sanity plot for gene: ", sanity_gene)

  p1 <- plot_exp(ref_sce, gene = sanity_gene, pt_size = 1.2) +
    ggtitle(paste0(sanity_gene, ": real data"))
  p2 <- plot_exp(sim_sce, gene = sanity_gene, pt_size = 1.2) +
    ggtitle(paste0(sanity_gene, ": simulated data"))
  save_plot(p1 + p2,
            file.path(sim_fig_dir, glue::glue("sanity_real_vs_sim_{sanity_gene}")),
            width = 10, height = 4)

  # ----------------------------------------------------------------------------
  # STEP 6: Generate shuffled (non-spatial) mean matrix
  # ----------------------------------------------------------------------------

  shuffle_idx <- sample(nrow(ref_para$mean_mat))
  non_de_mat <- ref_para$mean_mat[shuffle_idx, ]

  # ----------------------------------------------------------------------------
  # STEP 7: Alpha sweep - simulate at alpha = 0, 0.05, ..., 1.0
  # ----------------------------------------------------------------------------

  message("  Generating alpha sweep (21 levels)...")

  count <- lapply(seq(0, 1.0, 0.05), function(alpha) {
    sim_count <- simu_new(
      sce = ref_sce,
      mean_mat = alpha * ref_para$mean_mat + (1 - alpha) * non_de_mat,
      sigma_mat = ref_para$sigma_mat,
      zero_mat = ref_para$zero_mat,
      quantile_mat = NULL,
      copula_list = ref_copula$copula_list,
      n_cores = 5,
      family_use = "nb",
      input_data = ref_data$dat,
      new_covariate = ref_data$newCovariate,
      important_feature = rep(TRUE, dim(ref_sce)[1]),
      filtered_gene = NULL
    )
    rownames(sim_count) <- paste0(rownames(sim_count), "_", alpha)
    sim_count
  }) %>% do.call(rbind, .)

  # ----------------------------------------------------------------------------
  # STEP 8: Export counts and location matrices
  # ----------------------------------------------------------------------------

  write.csv(as.data.frame(ref_data$newCovariate),
            file = file.path(sim_data_dir, "location.csv"))
  write.csv(count, file = file.path(sim_data_dir, "counts.csv"))

  message("  Done: ", slice, " -> ", sim_data_dir)
}

message("\nAll slices processed.")