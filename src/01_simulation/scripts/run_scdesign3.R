library(Seurat)
library(SingleCellExperiment)
library(scDesign3)
library(scales)
library(ggplot2)
library(cowplot)
library(dplyr)

data_dir <- file.path("outputs", "scDesign3", "data")
figures_dir <- file.path("outputs", "scDesign3", "figures")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

ref_sce <- readRDS(file.path("data", "VISIUM_sce.rds"))

# Remove mitochondrial genes before fitting because they can dominate the
# spatial signal and add unnecessary noise to the downstream models.
mt_idx<- grep("mt-",rownames(ref_sce))
if(length(mt_idx)!=0){
    ref_sce <- ref_sce[-mt_idx,]
}

ref_sce

# Select a smaller set of spatially variable genes with Moran's I so the
# initial scDesign3 fitting step stays tractable on the reference dataset.
num_genes <- 100
loc <- colData(ref_sce)[, c("spatial1","spatial2")]
features <- FindSpatiallyVariableFeatures(counts(ref_sce), 
                                          spatial.location = loc, 
                                          selection.method = "moransi", 
                                          nfeatures = num_genes)

top.features <- features[order(features$p.value),]

top.features <- rownames(top.features[1:num_genes,])

de_idx <- which(rownames(ref_sce) %in% top.features)
ref_sce <- ref_sce[de_idx, ]

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


plot_exp <- function(sce, gene, pt_size=1){
  df_loc <- colData(sce)[, c("spatial1","spatial2")]
  df_exp <- as.data.frame(counts(sce)[gene, ])
  colnames(df_exp) <- c('exp')
  df_exp$exp <- rescale(log1p(df_exp$exp))
  
  df <- cbind(df_exp, df_loc)
  
  p <- ggplot(data = df, aes(x = .data$spatial1, y = .data$spatial2)) +
    geom_point(aes(x = .data$spatial1, y = .data$spatial2, color = .data$exp), size = pt_size) +
    scale_colour_gradientn(colors = viridis_pal(option = "magma")(10), limits=c(0, 1)) +
    theme_cowplot() +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    ggtitle(gene)
  
  return(p)
}

options(repr.plot.height = 4, repr.plot.width = 5)

for(gene in rownames(ref_sce)[1:4]){
  p <- plot_exp(ref_sce, gene = gene, pt_size = 1.2)
  print(p)
  save_plot(p, glue::glue("sanity_ref_{gene}"), width = 5, height = 4)
}

# Run the first scDesign3 pass on the reduced feature set. This step fits the
# marginal spatial model for each gene and is the main cost driver in the
# pipeline.
set.seed(2024)
# constructs the input data for fit_marginal.
ref_data <- construct_data(
  sce = ref_sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = c("spatial1", "spatial2"),
  other_covariates = NULL,
  corr_by = "1"
)

# fit expression of each gene with GP model
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

# Rank genes by deviance explained and keep the strongest spatial signals for
# the joint model. This trims the later copula and simulation steps.
dev_explain <- sapply(ref_marginal, function(x){
  sum = summary(x$fit)
  return(sum$dev.expl)
})
dev_ordered <- order(dev_explain, decreasing = TRUE)
num_de <- 50
ordered <- dev_explain[dev_ordered]
sel_genes <- names(ordered)[1:num_de]

ref_para <- extract_para(
  sce = ref_sce,
  marginal_list = ref_marginal,
  n_cores = 5,
  family_use = "nb",
  new_covariate = ref_data$newCovariate,
  data = ref_data$dat
)

dev_explain <- sapply(ref_marginal, function(x){
  sum = summary(x$fit)
  return(sum$dev.expl)
})
dev_ordered <- order(dev_explain, decreasing = TRUE)
num_de <- 50
ordered <- dev_explain[dev_ordered]
sel_genes <- names(ordered)[1:num_de]

ref_sce <- ref_sce[sel_genes, ]

options(repr.plot.height = 4, repr.plot.width = 10)

p1 <- plot_exp(ref_sce, gene='Ttr', pt_size = 1.2)+ 
  ggtitle("Ttr: real data")
p2 <- plot_exp(ref_sce, gene='Mbp', pt_size = 1.2)+ 
  ggtitle("Mbp: real data")

print(p1 + p2)
save_plot(p1 + p2, "sanity_ref_Ttr_Mbp", width = 10, height = 4)

# Rebuild the input object with the selected genes so the final joint model is
# fit only on the most informative features.
ref_data <- construct_data(
  sce = ref_sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = c("spatial1", "spatial2"),
  other_covariates = NULL,
  corr_by = "1"
)

# fit expression of each gene with GP model
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

# Fit the Gaussian copula to preserve gene-gene dependence structure across the
# selected genes before simulation.
ref_copula <- fit_copula(
  sce = ref_sce,
  assay_use = "counts",
  marginal_list = ref_marginal,
  family_use = "nb",
  copula = "gaussian",
  n_cores = 5,
  input_data = ref_data$dat
)

# Extract the fitted mean, variance, and zero-inflation parameters that will be
# reused to generate synthetic count matrices.
ref_para <- extract_para(
  sce = ref_sce,
  marginal_list = ref_marginal,
  n_cores = 5,
  family_use = "nb",
  new_covariate = ref_data$newCovariate,
  data = ref_data$dat
)

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
  filtered_gene = NULL)

sim_sce <- SingleCellExperiment(list(counts = sim_count), 
                                colData = ref_data$newCovariate)

options(repr.plot.height = 4, repr.plot.width = 10)

p1 <- plot_exp(ref_sce, gene='Ttr', pt_size = 1.2)+ 
  ggtitle("Ttr: real data")
p2 <- plot_exp(sim_sce, gene='Ttr', pt_size = 1.2)+ 
  ggtitle("Ttr: simulated data")

print(p1 + p2)
save_plot(p1 + p2, "sanity_real_vs_sim_Ttr", width = 10, height = 4)

# Generate a shuffled mean matrix to break the spatial pattern while keeping the
# marginal distribution shape available for controlled mixing experiments.
shuffle_idx <- sample(nrow(ref_para$mean_mat))
non_de_mat <- ref_para$mean_mat[shuffle_idx, ]

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
  filtered_gene = NULL)

sim_sce1 <- SingleCellExperiment(list(counts =sim_count1), 
                                 colData = ref_data$newCovariate)

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
  filtered_gene = NULL)

sim_sce2 <- SingleCellExperiment(list(counts =sim_count2), 
                                 colData = ref_data$newCovariate)

options(repr.plot.height = 8, repr.plot.width = 10)

p1 <- plot_exp(ref_sce, gene='Nrgn', pt_size = 1.1) + 
  ggtitle("Nrgn: real data")
p2 <- plot_exp(sim_sce, gene='Nrgn', pt_size = 1.1) + 
  ggtitle("Nrgn: simualte data 100% signal")
p3 <- plot_exp(sim_sce1, gene='Nrgn', pt_size = 1.1) + 
  ggtitle("Nrgn: simualte data 90% signal")
p4 <- plot_exp(sim_sce2, gene='Nrgn', pt_size = 1.1) + 
  ggtitle("Nrgn: simualte data 10% signal")

print(p1 + p2 + p3 + p4)
save_plot(p1 + p2 + p3 + p4, "compare_Nrgn", width = 10, height = 8)

options(repr.plot.height = 8, repr.plot.width = 10)

p1 <- plot_exp(ref_sce, gene='Ttr', pt_size = 1.1) + 
  ggtitle("Ttr: real data")
p2 <- plot_exp(sim_sce, gene='Ttr', pt_size = 1.1) + 
  ggtitle("Ttr: simualte data 100% signal")
p3 <- plot_exp(sim_sce1, gene='Ttr', pt_size = 1.1) + 
  ggtitle("Ttr: simualte data 90% signal")
p4 <- plot_exp(sim_sce2, gene='Ttr', pt_size = 1.1) + 
  ggtitle("Ttr: simualte data 10% signal")

print(p1 + p2 + p3 + p4)
save_plot(p1 + p2 + p3 + p4, "compare_Ttr", width = 10, height = 8)

options(repr.plot.height = 8, repr.plot.width = 10)

p1 <- plot_exp(ref_sce, gene='S100a5', pt_size = 1.1) + 
  ggtitle("S100a5: real data")
p2 <- plot_exp(sim_sce, gene='S100a5', pt_size = 1.1) + 
  ggtitle("S100a5: simualte data 100% signal")
p3 <- plot_exp(sim_sce1, gene='S100a5', pt_size = 1.1) + 
  ggtitle("S100a5: simualte data 90% signal")
p4 <- plot_exp(sim_sce2, gene='S100a5', pt_size = 1.1) + 
  ggtitle("S100a5: simualte data 10% signal")

print(p1 + p2 + p3 + p4)
save_plot(p1 + p2 + p3 + p4, "compare_S100a5", width = 10, height = 8)

options(repr.plot.height = 8, repr.plot.width = 10)

p1 <- plot_exp(ref_sce, gene='Doc2g', pt_size = 1.1) + 
  ggtitle("Doc2g: real data")
p2 <- plot_exp(sim_sce, gene='Doc2g', pt_size = 1.1) + 
  ggtitle("Doc2g: simualte data 100% signal")
p3 <- plot_exp(sim_sce1, gene='Doc2g', pt_size = 1.1) + 
  ggtitle("Doc2g: simualte data 90% signal")
p4 <- plot_exp(sim_sce2, gene='Doc2g', pt_size = 1.1) + 
  ggtitle("Doc2g: simualte data 10% signal")

print(p1 + p2 + p3 + p4)
save_plot(p1 + p2 + p3 + p4, "compare_Doc2g", width = 10, height = 8)

# Sweep across mixing weights to generate datasets with progressively weaker
# spatial signal and quantify how the benchmark behaves across the full range.
count <- lapply(seq(0, 1.0, 0.2), function(alpha){
  sim_count <- simu_new(
    sce = ref_sce,
    mean_mat =  alpha * ref_para$mean_mat + (1 - alpha) * non_de_mat,
    sigma_mat = ref_para$sigma_mat,
    zero_mat = ref_para$zero_mat,
    quantile_mat = NULL,
    copula_list = ref_copula$copula_list,
    n_cores = 1,
    family_use = "nb",
    input_data = ref_data$dat,
    new_covariate = ref_data$newCovariate,
    important_feature = rep(TRUE, dim(ref_sce)[1]),
    filtered_gene = NULL)
  
  rownames(sim_count) <- paste0(rownames(sim_count), "_", alpha)
  
  return(sim_count)
  
}) %>% do.call(rbind, .)

# Save the final benchmark artifacts: coordinates and the stacked simulated
# count matrix across all alpha settings.
saveRDS(sim_sce, file = glue::glue("{data_dir}/sim_sce.rds"))
write.csv(ref_data$newCovariate, file = glue::glue("{data_dir}/location.csv"))
write.csv(count, file = glue::glue('{data_dir}/counts.csv'))