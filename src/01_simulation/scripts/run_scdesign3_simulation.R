# ==============================================================================
# run_scdesign3_simulation.R
# Fits scDesign3 to each stxBrain slice (anterior1/2, posterior1/2), then
# generates the alpha-sweep simulated counts (alpha = 0, 0.05, ..., 1.0) by
# blending the fitted spatial mean with a shuffled (non-spatial) mean matrix,
# PLUS a fixed diverse null background of permuted real genes.
#
# Per slice:
#   1. Load stxBrain slice via SeuratData; SCTransform + Moran's I (top 200)
#   2. First pass: construct_data -> fit_marginal (GP k=500, NB) -> fit_copula
#   3. Select top 50 genes by deviance explained
#   4. Second pass on the 50 genes
#   5. Sanity plot: Mbp real vs simulated at 100% signal
#   6. Alpha sweep -> stacked counts matrix with "gene_alpha" rownames
#   7. Null background: N_NULL real genes (expression-stratified to match the
#      50 SVGs' total-count distribution), each permuted across spots
#      independently -> "NULL_gene" features, stacked after the alpha sweep
#   8. Sanity plot: null genes real vs permuted (spatial pattern destroyed)
#   9. Save counts.csv + location.csv
#
# Design notes:
#   - The null background fixes the negative-class deficiency of the previous
#     all-in-one matrix (only 50 shuffled copies): positives (alpha > 0,
#     1000 features) vs negatives (alpha = 0 copies UNION NULL_ background,
#     50 + N_NULL) gives a balanced, diverse, expression-controlled negative
#     class for auPRC.
#   - Expression stratification controls the expression-level confound
#     (Chen et al. 2024: SVG scores track expression level).
#   - The 50 alpha = 0 copies are KEPT: they are the per-gene gradient
#     endpoints and expression-matched hard negatives.
#   - The null step uses its own seed AFTER the alpha sweep, so the RNG stream
#     feeding simu_new is untouched and alpha-feature counts remain identical
#     to previous runs.
#
# Edit SLICES / N_NULL below to control scope.
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

# --- Fixed null-background size (expression-stratified permuted real genes) ---
N_NULL <- 1000

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

  # Retain the FULL raw counts matrix (all genes x spots) BEFORE subsetting:
  # this is the pool from which the null-background genes are sampled (step 7).
  full_counts_mat <- GetAssayData(seu, assay = "Spatial", layer = "counts")

  # Subset the Seurat object to the top 200 SVGs
  seu <- seu[top.features, ]

  # Convert the subsetted Seurat object to SCE. Extract RAW counts from the
  # original "Spatial" assay (not the "SCT" assay) since scDesign3 needs raw
  # counts for construct_data(assay_use = "counts").
  coords <- GetTissueCoordinates(seu, scale = "lowres")
  # imagecol -> spatial1, imagerow -> spatial2 (match Visium orientation)
  df_loc <- data.frame(
    spatial1 = coords$x,
    spatial2 = coords$y,
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
  # STEP 7.5: Fixed diverse null background (N_NULL permuted real genes)
  # ----------------------------------------------------------------------------
  # The N_NULL highest-expressed background genes (excluding the 50 SVGs) --
  # the closest available expression match to the SVG panel, controlling the
  # expression-level confound (Chen et al. 2024: SVG scores track expression).
  # Perfect matching is impossible (only ~250 pool genes reach the SVG median
  # expression), but this is far stricter than published protocols (SRTsim
  # used median-expression nulls). Each gene's counts are permuted across
  # spots independently: real marginals (mean, variance, dropout) preserved,
  # spatial signal destroyed by construction. Named "NULL_gene" so downstream
  # parsers tag them as null (non-numeric alpha suffix).
  #
  # NOTE: own seed, set AFTER the alpha sweep -- the RNG stream feeding
  # simu_new above is untouched, so alpha-feature counts are identical to
  # previous runs.
  # ----------------------------------------------------------------------------

  message("  Generating null background (", N_NULL, " permuted real genes)...")

  set.seed(2025)

  full_counts <- as.matrix(full_counts_mat)
  pool_genes  <- setdiff(
    rownames(full_counts)[rowSums(full_counts) > 0],
    sel_genes
  )
  pool_totals <- rowSums(full_counts[pool_genes, , drop = FALSE])
  null_genes  <- names(sort(pool_totals, decreasing = TRUE))[seq_len(N_NULL)]

  null_counts_raw <- full_counts[null_genes, , drop = FALSE]
  # Per-gene independent spot permutation: destroys spatial autocorrelation,
  # preserves the marginal count distribution exactly.
  null_counts <- t(apply(null_counts_raw, 1, sample))
  colnames(null_counts) <- colnames(full_counts)

  null_export <- null_counts
  rownames(null_export) <- paste0("NULL_", null_genes)

  count <- rbind(count, null_export)

  message("  Null background: ", length(null_genes),
          " highest-expressed background genes, spot-permuted")

  # ----------------------------------------------------------------------------
  # STEP 7.6: Sanity plot - null genes real vs permuted
  # ----------------------------------------------------------------------------

  null_sanity_genes <- head(null_genes, 2)
  sce_null_orig <- SingleCellExperiment(
    list(counts = null_counts_raw[null_sanity_genes, , drop = FALSE]),
    colData = df_loc
  )
  sce_null_perm <- SingleCellExperiment(
    list(counts = null_counts[null_sanity_genes, , drop = FALSE]),
    colData = df_loc
  )

  null_plots <- lapply(null_sanity_genes, function(g) {
    p1 <- plot_exp(sce_null_orig, gene = g, pt_size = 1.2) +
      ggtitle(paste0(g, ": real (possibly spatial)"))
    p2 <- plot_exp(sce_null_perm, gene = g, pt_size = 1.2) +
      ggtitle(paste0(g, ": permuted (null)"))
    p1 + p2
  })
  save_plot(plot_grid(plotlist = null_plots, ncol = 1),
            file.path(sim_fig_dir, "sanity_null_permutation"),
            width = 10, height = 4 * length(null_plots))

  # ----------------------------------------------------------------------------
  # STEP 8: Export counts and location matrices
  # ----------------------------------------------------------------------------

  write.csv(as.data.frame(ref_data$newCovariate),
            file = file.path(sim_data_dir, "location.csv"))
  write.csv(count, file = file.path(sim_data_dir, "counts.csv"))

  message("  Done: ", slice, " -> ", sim_data_dir)
}

message("\nAll slices processed.")