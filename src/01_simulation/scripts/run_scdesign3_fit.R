# ==============================================================================
# run_scdesign3_fit.R
# Fits the scDesign3 generative model to a reference Visium dataset, then
# generates simulated counts at multiple signal-strength levels (100%, 90%, 10%)
# for downstream SVG rotation invariance benchmarking.
#
# Pipeline overview:
#   1. Load and preprocess reference Visium SCE
#   2. Select top spatially variable features (Moran's I)
#   3. Fit scDesign3: marginal distributions -> copula -> parameter extraction
#   4. Simulate full-signal data and validate against real data
#   5. Create non-spatial (shuffled) mean matrix for signal dilution
#   6. Simulate at 90% and 10% signal levels for validation
#   7. Export fitted model for use by run_scdesign3_sweep.R
#
# Output: scdesign3_fitted.rds (ref_sce, ref_para, ref_copula, ref_data, non_de_mat)
# ==============================================================================

library(Seurat)
library(SingleCellExperiment)
library(scDesign3)
library(scales)
library(ggplot2)
library(cowplot)
library(dplyr)

# --- Setup: resolve project root and output directories ---
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

data_dir <- file.path(project_root, "src", "01_simulation", "outputs", "scDesign3", "data")
figures_dir <- file.path(project_root, "src", "01_simulation", "outputs", "scDesign3", "figures")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# SECTION 1: LOAD AND PREPROCESS REFERENCE DATA
# ==============================================================================

# Load reference Visium SCE (containing counts + spatial coordinates)
ref_sce <- readRDS(file.path(project_root, "data", "VISIUM_sce.rds"))

# Remove mitochondrial genes (prefixed "mt-") to avoid confounding
mt_idx <- grep("mt-", rownames(ref_sce))
if (length(mt_idx) != 0) {
  ref_sce <- ref_sce[-mt_idx, ]
}

ref_sce

# ==============================================================================
# SECTION 2: SELECT TOP SPATIALLY VARIABLE FEATURES
# ==============================================================================

# Use Moran's I (via Seurat) to identify the top 200 SVGs in the reference data.
# These genes form the basis for the simulation.
num_genes <- 200
loc <- colData(ref_sce)[, c("spatial1", "spatial2")]
features <- FindSpatiallyVariableFeatures(counts(ref_sce),
  spatial.location = loc,
  selection.method = "moransi",
  nfeatures = num_genes
)

# Order by significance and keep the top 200
top.features <- features[order(features$p.value), ]
top.features <- rownames(top.features[1:num_genes, ])

# Subset the SCE to only these SVGs
de_idx <- which(rownames(ref_sce) %in% top.features)
ref_sce <- ref_sce[de_idx, ]

# ==============================================================================
# SECTION 3: SANITY PLOT HELPERS
# ==============================================================================

# Helper to save plots as both PNG and PDF
save_plot <- function(plot_object, file_stub, width = 6, height = 5) {
  ggsave(
    filename = glue::glue("{figures_dir}/{file_stub}.png"),
    plot = plot_object,
    width = width,
    height = height,
    dpi = 300
  )
  ggsave(
    filename = glue::glue("{figures_dir}/{file_stub}.pdf"),
    plot = plot_object,
    width = width,
    height = height
  )
}

# Helper to plot spatial expression of a single gene
plot_exp <- function(sce, gene, pt_size = 1) {
  df_loc <- colData(sce)[, c("spatial1", "spatial2")]
  df_exp <- as.data.frame(counts(sce)[gene, ])
  colnames(df_exp) <- c("exp")
  # Rescale log1p expression to [0, 1] for consistent color scale
  df_exp$exp <- rescale(log1p(df_exp$exp))

  df <- cbind(df_exp, df_loc)

  p <- ggplot(data = df, aes(x = .data$spatial1, y = .data$spatial2)) +
    geom_point(aes(x = .data$spatial1, y = .data$spatial2, color = .data$exp), size = pt_size) +
    scale_colour_gradientn(colors = viridis_pal(option = "magma")(10), limits = c(0, 1)) +
    theme_cowplot() +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    ggtitle(gene)

  return(p)
}

options(repr.plot.height = 4, repr.plot.width = 5)

# --- Sanity check: plot first 4 reference genes ---
for (gene in rownames(ref_sce)[1:4]) {
  p <- plot_exp(ref_sce, gene = gene, pt_size = 1.2)
  print(p)
  save_plot(p, glue::glue("sanity_ref_{gene}"), width = 5, height = 4)
}

# ==============================================================================
# SECTION 4: scDesign3 MODEL FITTING (FIRST PASS ON FULL GENE SET)
# ==============================================================================

set.seed(2024)

# Step 1: Construct scDesign3 data object from SCE
#   celltype = "cell_type": account for cell type as a covariate
#   spatial = c("spatial1", "spatial2"): model spatial coordinates
#   corr_by = "1": assume all spots share one correlation group
ref_data <- construct_data(
  sce = ref_sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = c("spatial1", "spatial2"),
  other_covariates = NULL,
  corr_by = "1"
)

# Step 2: Fit marginal distributions (negative binomial GAM with GP smooth)
#   mu_formula: GP smooth over spatial coordinates (k = 500 basis functions)
#   sigma_formula: constant dispersion
#   family_use = "nb": negative binomial
ref_marginal <- fit_marginal(
  data = ref_data,
  predictor = "gene",
  mu_formula = "s(spatial1, spatial2, bs = 'gp', k = 500)",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 5,
  usebam = FALSE,
  trace = TRUE
)

# Step 3: Select genes with highest deviance explained (best spatial fit)
# We keep the top 50 genes for the final simulation to reduce dimensionality
dev_explain <- sapply(ref_marginal, function(x) {
  sum <- summary(x$fit)
  return(sum$dev.expl)
})
dev_ordered <- order(dev_explain, decreasing = TRUE)
num_de <- 50
ordered <- dev_explain[dev_ordered]
sel_genes <- names(ordered)[1:num_de]

# Step 4: Extract marginal parameters (mean, sigma, zero-inflation matrices)
ref_para <- extract_para(
  sce = ref_sce,
  marginal_list = ref_marginal,
  n_cores = 5,
  family_use = "nb",
  new_covariate = ref_data$newCovariate,
  data = ref_data$dat
)

# --- Re-derive deviance and re-select top 50 genes for consistency ---
dev_explain <- sapply(ref_marginal, function(x) {
  sum <- summary(x$fit)
  return(sum$dev.expl)
})
dev_ordered <- order(dev_explain, decreasing = TRUE)
num_de <- 50
ordered <- dev_explain[dev_ordered]
sel_genes <- names(ordered)[1:num_de]

# Subset reference SCE to the 50 best-fit genes
ref_sce <- ref_sce[sel_genes, ]

options(repr.plot.height = 4, repr.plot.width = 10)

# --- Sanity check: Ttr and Mbp spatial expression in real data ---
p1 <- plot_exp(ref_sce, gene = "Ttr", pt_size = 1.2) +
  ggtitle("Ttr: real data")
p2 <- plot_exp(ref_sce, gene = "Mbp", pt_size = 1.2) +
  ggtitle("Mbp: real data")

print(p1 + p2)
save_plot(p1 + p2, "sanity_ref_Ttr_Mbp", width = 10, height = 4)

# ==============================================================================
# SECTION 5: scDesign3 MODEL FITTING (SECOND PASS ON TOP-50 GENES)
# ==============================================================================

# Reconstruct data object with the reduced gene set
ref_data <- construct_data(
  sce = ref_sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = c("spatial1", "spatial2"),
  other_covariates = NULL,
  corr_by = "1"
)

# Re-fit marginal distributions on the reduced gene set
ref_marginal <- fit_marginal(
  data = ref_data,
  predictor = "gene",
  mu_formula = "s(spatial1, spatial2, bs = 'gp', k = 500)",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 5,
  usebam = FALSE,
  trace = TRUE
)

# Step 5: Fit Gaussian copula to capture gene-gene correlations
ref_copula <- fit_copula(
  sce = ref_sce,
  assay_use = "counts",
  marginal_list = ref_marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 5,
  input_data = ref_data$dat
)

# Step 6: Extract final parameters (mean, sigma, zero-inflation matrices)
ref_para <- extract_para(
  sce = ref_sce,
  marginal_list = ref_marginal,
  n_cores = 5,
  family_use = "nb",
  new_covariate = ref_data$newCovariate,
  data = ref_data$dat
)

# ==============================================================================
# SECTION 6: SIMULATE AT 100% SIGNAL AND VALIDATE
# ==============================================================================

# Generate synthetic counts at 100% signal strength (alpha = 1.0)
sim_count <- simu_new(
  sce = ref_sce,
  mean_mat = ref_para$mean_mat,
  sigma_mat = ref_para$sigma_mat,
  zero_mat = ref_para$zero_mat,
  quantile_mat = NULL,
  copula_list = ref_copula$copula_list,
  n_cores = 1,
  family_use = "nb",
  input_data = ref_data$dat,
  new_covariate = ref_data$newCovariate,
  important_feature = rep(TRUE, dim(ref_sce)[1]),
  filtered_gene = NULL
)

# Wrap simulated counts in a SingleCellExperiment
sim_sce <- SingleCellExperiment(list(counts = sim_count),
  colData = ref_data$newCovariate
)

options(repr.plot.height = 4, repr.plot.width = 10)

# --- Validate: real vs simulated Ttr ---
p1 <- plot_exp(ref_sce, gene = "Ttr", pt_size = 1.2) +
  ggtitle("Ttr: real data")
p2 <- plot_exp(sim_sce, gene = "Ttr", pt_size = 1.2) +
  ggtitle("Ttr: simulated data")

print(p1 + p2)
save_plot(p1 + p2, "sanity_real_vs_sim_Ttr", width = 10, height = 4)

# ==============================================================================
# SECTION 7: CREATE NON-SPATIAL MEAN MATRIX (SHUFFLED BACKGROUND)
# ==============================================================================

# Shuffle the mean matrix rows to break spatial structure.
# This creates a non-DE baseline for signal dilution in the alpha sweep.
shuffle_idx <- sample(nrow(ref_para$mean_mat))
non_de_mat <- ref_para$mean_mat[shuffle_idx, ]

# ==============================================================================
# SECTION 8: SIMULATE AT 90% AND 10% SIGNAL FOR VALIDATION
# ==============================================================================

# --- 90% signal (alpha = 0.9), 10% non-spatial background ---
sim_count1 <- simu_new(
  sce = ref_sce,
  mean_mat = 0.9 * ref_para$mean_mat + 0.1 * non_de_mat,
  sigma_mat = ref_para$sigma_mat,
  zero_mat = ref_para$zero_mat,
  quantile_mat = NULL,
  copula_list = ref_copula$copula_list,
  n_cores = 1,
  family_use = "nb",
  input_data = ref_data$dat,
  new_covariate = ref_data$newCovariate,
  important_feature = rep(TRUE, dim(ref_sce)[1]),
  filtered_gene = NULL
)

sim_sce1 <- SingleCellExperiment(list(counts = sim_count1),
  colData = ref_data$newCovariate
)

# --- 10% signal (alpha = 0.1), 90% non-spatial background ---
sim_count2 <- simu_new(
  sce = ref_sce,
  mean_mat = 0.1 * ref_para$mean_mat + 0.9 * non_de_mat,
  sigma_mat = ref_para$sigma_mat,
  zero_mat = ref_para$zero_mat,
  quantile_mat = NULL,
  copula_list = ref_copula$copula_list,
  n_cores = 1,
  family_use = "nb",
  input_data = ref_data$dat,
  new_covariate = ref_data$newCovariate,
  important_feature = rep(TRUE, dim(ref_sce)[1]),
  filtered_gene = NULL
)

sim_sce2 <- SingleCellExperiment(list(counts = sim_count2),
  colData = ref_data$newCovariate
)

# ==============================================================================
# SECTION 9: COMPARE ACROSS SIGNAL LEVELS FOR 4 EXAMPLE GENES
# ==============================================================================

options(repr.plot.height = 8, repr.plot.width = 10)

# --- Nrgn: real vs 100% vs 90% vs 10% ---
p1 <- plot_exp(ref_sce, gene = "Nrgn", pt_size = 1.1) +
  ggtitle("Nrgn: real data")
p2 <- plot_exp(sim_sce, gene = "Nrgn", pt_size = 1.1) +
  ggtitle("Nrgn: simulated data 100% signal")
p3 <- plot_exp(sim_sce1, gene = "Nrgn", pt_size = 1.1) +
  ggtitle("Nrgn: simulated data 90% signal")
p4 <- plot_exp(sim_sce2, gene = "Nrgn", pt_size = 1.1) +
  ggtitle("Nrgn: simulated data 10% signal")

print(p1 + p2 + p3 + p4)
save_plot(p1 + p2 + p3 + p4, "compare_Nrgn", width = 10, height = 8)

options(repr.plot.height = 8, repr.plot.width = 10)

# --- Ttr: real vs 100% vs 90% vs 10% ---
p1 <- plot_exp(ref_sce, gene = "Ttr", pt_size = 1.1) +
  ggtitle("Ttr: real data")
p2 <- plot_exp(sim_sce, gene = "Ttr", pt_size = 1.1) +
  ggtitle("Ttr: simulated data 100% signal")
p3 <- plot_exp(sim_sce1, gene = "Ttr", pt_size = 1.1) +
  ggtitle("Ttr: simulated data 90% signal")
p4 <- plot_exp(sim_sce2, gene = "Ttr", pt_size = 1.1) +
  ggtitle("Ttr: simulated data 10% signal")

print(p1 + p2 + p3 + p4)
save_plot(p1 + p2 + p3 + p4, "compare_Ttr", width = 10, height = 8)

options(repr.plot.height = 8, repr.plot.width = 10)

# --- S100a5: real vs 100% vs 90% vs 10% ---
p1 <- plot_exp(ref_sce, gene = "S100a5", pt_size = 1.1) +
  ggtitle("S100a5: real data")
p2 <- plot_exp(sim_sce, gene = "S100a5", pt_size = 1.1) +
  ggtitle("S100a5: simulated data 100% signal")
p3 <- plot_exp(sim_sce1, gene = "S100a5", pt_size = 1.1) +
  ggtitle("S100a5: simulated data 90% signal")
p4 <- plot_exp(sim_sce2, gene = "S100a5", pt_size = 1.1) +
  ggtitle("S100a5: simulated data 10% signal")

print(p1 + p2 + p3 + p4)
save_plot(p1 + p2 + p3 + p4, "compare_S100a5", width = 10, height = 8)

options(repr.plot.height = 8, repr.plot.width = 10)

# --- Doc2g: real vs 100% vs 90% vs 10% ---
p1 <- plot_exp(ref_sce, gene = "Doc2g", pt_size = 1.1) +
  ggtitle("Doc2g: real data")
p2 <- plot_exp(sim_sce, gene = "Doc2g", pt_size = 1.1) +
  ggtitle("Doc2g: simulated data 100% signal")
p3 <- plot_exp(sim_sce1, gene = "Doc2g", pt_size = 1.1) +
  ggtitle("Doc2g: simulated data 90% signal")
p4 <- plot_exp(sim_sce2, gene = "Doc2g", pt_size = 1.1) +
  ggtitle("Doc2g: simulated data 10% signal")

print(p1 + p2 + p3 + p4)
save_plot(p1 + p2 + p3 + p4, "compare_Doc2g", width = 10, height = 8)

# ==============================================================================
# SECTION 10: EXPORT FITTED MODEL
# ==============================================================================

# Save all model components needed by run_scdesign3_sweep.R for the alpha sweep
saveRDS(
  list(
    ref_sce    = ref_sce,      # Reference SCE (top-50 SVGs)
    ref_para   = ref_para,     # Marginal parameters (mean, sigma, zero)
    ref_copula = ref_copula,   # Gaussian copula (gene-gene correlations)
    ref_data   = ref_data,     # scDesign3 data object
    non_de_mat = non_de_mat    # Shuffled mean matrix for signal dilution
  ),
  file = file.path(data_dir, "scdesign3_fitted.rds")
)
