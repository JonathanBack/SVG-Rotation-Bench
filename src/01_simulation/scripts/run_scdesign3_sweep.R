# ==============================================================================
# run_scdesign3_sweep.R
# Generates simulated spatial transcriptomics count data across a sweep of
# signal-strength values (alpha = 0, 0.05, ..., 1.0) by blending the fitted
# scDesign3 mean matrix with a shuffled (non-spatial) mean matrix.
# Output: location.csv and counts.csv consumed by downstream rotation/benchmark.
# ==============================================================================

library(SingleCellExperiment)
library(scDesign3)
library(dplyr)

# --- Setup: resolve project root and output directory ---
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

data_dir <- file.path(project_root, "src", "01_simulation", "outputs", "scDesign3", "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# Load the pre-fitted scDesign3 model (ref_sce, ref_para, ref_copula, ref_data, non_de_mat)
fitted <- readRDS(file.path(data_dir, "scdesign3_fitted.rds"))
list2env(fitted, .GlobalEnv)

# --- Alpha sweep: blend true spatial signal with shuffled (non-spatial) background
# For each alpha in 0, 0.05, ..., 1.0, generate count using:
#   mean_mat = alpha * ref_para$mean_mat + (1 - alpha) * non_de_mat
# Features are named "gene_alpha" so the alpha value is encoded in the rowname.
count <- lapply(seq(0, 1.0, 0.05), function(alpha) {
  sim_count <- simu_new(
    sce = ref_sce,
    mean_mat = alpha * ref_para$mean_mat + (1 - alpha) * non_de_mat,
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

  rownames(sim_count) <- paste0(rownames(sim_count), "_", alpha)
  sim_count
}) %>% do.call(rbind, .)

# --- Export: spatial coordinates and stacked count matrix ---
write.csv(as.data.frame(ref_data$newCovariate), file = file.path(data_dir, "location.csv"))
write.csv(count, file = file.path(data_dir, "counts.csv"))
