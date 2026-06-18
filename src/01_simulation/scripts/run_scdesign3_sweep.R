library(SingleCellExperiment)
library(scDesign3)
library(dplyr)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

data_dir <- file.path(project_root, "src", "01_simulation", "outputs", "scDesign3", "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

fitted <- readRDS(file.path(data_dir, "scdesign3_fitted.rds"))
list2env(fitted, .GlobalEnv)

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

write.csv(as.data.frame(ref_data$newCovariate), file = file.path(data_dir, "location.csv"))
write.csv(count, file = file.path(data_dir, "counts.csv"))
