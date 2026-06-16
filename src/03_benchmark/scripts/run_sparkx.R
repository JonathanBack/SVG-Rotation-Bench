library(SPARK)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

counts_file <- file.path(project_root, "src", "01_simulation", "outputs", "scDesign3", "data", "counts.csv")
locations_template <- file.path(project_root, "src", "02_rotation", "outputs", "locations", "rotated_locations_%s.csv")
output_dir <- file.path(project_root, "src", "03_benchmark", "outputs", "sparkx")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
counts <- counts[!duplicated(rownames(counts)), ]

angles_degrees <- c(0, 30, 45, 60)

for (angle in angles_degrees) {
  output_file <- file.path(output_dir, paste0("scdesign3_angle", angle, ".csv"))

  if (file.exists(output_file)) {
    message("Skipping angle ", angle, " -- output already exists: ", output_file)
    next
  }

  message("Running SPARK-X for angle = ", angle)

  location_file <- sprintf(locations_template, angle)
  locations <- read.csv(location_file, row.names = 1, check.names = FALSE)
  locs <- as.matrix(locations[colnames(counts), c("spatial1", "spatial2")])

  result <- sparkx(counts, locs, numCores = 1, option = "mixture")

  df_out <- data.frame(
    gene = rownames(result$res_mtest),
    combined_pval = result$res_mtest[, 1],
    adjusted_pval_BY = result$res_mtest[, 2],
    row.names = NULL,
    check.names = FALSE
  )
  write.csv(df_out, output_file, row.names = FALSE)

  message("Saved ", output_file)
}

message("SPARK-X benchmark complete.")
