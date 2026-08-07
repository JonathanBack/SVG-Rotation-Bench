# ==============================================================================
# sanity_simulation_plots.R
# Standalone visualization of simulated spatial expression across alpha levels
# for a chosen gene. Plots the gene at alpha = 1, 0.6, 0.2, 0 side-by-side so
# that the dilution of spatial signal can be visually verified.
#
# Edit SLICES and GENE below to control scope and target gene.
# Output: src/01_simulation/outputs/scDesign3/{slice}/simulated/figures/
# ==============================================================================

library(SingleCellExperiment)
library(scales)
library(ggplot2)
library(cowplot)
library(patchwork)

# --- Select which slices and which gene to visualize ---
SLICES <- c("anterior1", "anterior2", "posterior1", "posterior2")
GENE   <- "Camk2n1"  # if absent for a slice, falls back to first available

ALPHA_LEVELS <- c("1", "0.9", "0.6", "0.3", "0")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

plot_exp <- function(loc_df, exp_vec, gene, alpha_label, pt_size = 1.2) {
  df_exp <- data.frame(exp = exp_vec)
  df_exp$exp <- rescale(log1p(df_exp$exp))
  df <- cbind(df_exp, loc_df)

  ggplot(data = df, aes(x = .data$spatial1, y = .data$spatial2)) +
    geom_point(aes(color = .data$exp), size = pt_size) +
    scale_colour_gradientn(colors = viridis_pal(option = "magma")(10),
                           limits = c(0, 1)) +
    theme_cowplot() +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    ggtitle(paste0(gene, " | alpha = ", alpha_label))
}

for (slice in SLICES) {

  sim_data_dir <- file.path(project_root, "src", "01_simulation", "outputs",
                            "scDesign3", slice, "simulated", "data")
  sim_fig_dir  <- file.path(project_root, "src", "01_simulation", "outputs",
                            "scDesign3", slice, "simulated", "figures")
  dir.create(sim_fig_dir, recursive = TRUE, showWarnings = FALSE)

  counts_file <- file.path(sim_data_dir, "counts.csv")
  loc_file    <- file.path(sim_data_dir, "location.csv")
  if (!file.exists(counts_file) || !file.exists(loc_file)) {
    message("Skipping ", slice, " -- simulated data not found")
    next
  }
  message("Building sanity simulation plots for slice: ", slice)

  counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
  loc_df <- read.csv(loc_file, row.names = 1, check.names = FALSE)
  loc_df <- loc_df[, c("spatial1", "spatial2")]

  # Determine target gene (falls back if GENE absent for this slice)
  available <- unique(sapply(strsplit(rownames(counts), "_"), `[`, 1))
  target_gene <- if (GENE %in% available) GENE else available[1]
  if (!GENE %in% available)
    message("  Gene ", GENE, " not present in slice ", slice,
            " -- falling back to ", target_gene)

  # Pick the rows matching the target gene at each alpha level
  plots <- list()
  for (alpha_label in ALPHA_LEVELS) {
    feature_name <- paste0(target_gene, "_", alpha_label)
    if (!feature_name %in% rownames(counts)) {
      message("  Missing feature: ", feature_name, " -- skipping panel")
      next
    }
    plots[[alpha_label]] <- plot_exp(
      loc_df, as.numeric(counts[feature_name, ]),
      gene = target_gene, alpha_label = alpha_label
    )
  }

  if (length(plots) == 0) {
    message("  No matching features found for ", target_gene, " in slice ", slice)
    next
  }

  combined <- wrap_plots(plots, ncol = length(plots))
  ggsave(
    filename = file.path(sim_fig_dir,
                         glue::glue("sanity_simulation_{target_gene}.png")),
    plot = combined, width = 4 * length(plots), height = 4, dpi = 300
  )
  ggsave(
    filename = file.path(sim_fig_dir,
                         glue::glue("sanity_simulation_{target_gene}.pdf")),
    plot = combined, width = 4 * length(plots), height = 4
  )
}

message("\nSanity simulation plots complete.")