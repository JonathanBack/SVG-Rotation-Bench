# ==============================================================================
# run_sparkx.R
# Benchmarks the SPARK-X method for SVG detection across rotated datasets.
# For each rotation angle, loads the simulated counts (shared across angles) and
# the rotated spatial locations, runs SPARK-X, and saves results.
# Output: scdesign3_angle{angle}_results.rds and runtime CSV.
# ==============================================================================

library(SPARK)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# --- Setup: paths to shared counts and per-angle rotated locations ---
counts_file <- file.path(project_root, "src", "01_simulation", "outputs", "scDesign3", "data", "counts.csv")
locations_template <- file.path(project_root, "src", "02_rotation", "outputs", "locations", "rotated_locations_%s.csv")
output_dir <- file.path(project_root, "src", "03_benchmark", "outputs", "sparkx")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load count matrix (genes x cells); remove duplicate gene names if any
counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
counts <- counts[!duplicated(rownames(counts)), ]

angles_degrees <- c(0, 30, 45, 60)

for (angle in angles_degrees) {
  rds_file <- file.path(output_dir, paste0("scdesign3_angle", angle, "_results.rds"))
  runtime_file <- file.path(output_dir, paste0("scdesign3_angle", angle, "_runtime.csv"))

  if (file.exists(rds_file) && file.exists(runtime_file)) {
    message("Skipping angle ", angle, " -- outputs already exist")
    next
  }

  message("Running SPARK-X for angle = ", angle)

  # Load the rotated spatial coordinates for this angle
  location_file <- sprintf(locations_template, angle)
  locations <- read.csv(location_file, row.names = 1, check.names = FALSE)
  # Subset and order locations to match counts column order
  locs <- as.matrix(locations[colnames(counts), c("spatial1", "spatial2")])

  # --- Run SPARK-X with mixture model option ---
  t_start <- proc.time()
  result <- sparkx(counts, locs, numCores = 1, option = "mixture")
  elapsed <- proc.time() - t_start

  saveRDS(result, rds_file)
  write.csv(
    data.frame(angle = angle, elapsed_sec = unname(elapsed["elapsed"])),
    runtime_file,
    row.names = FALSE
  )

  message("Saved ", rds_file, " (", round(elapsed["elapsed"], 1), "s)")
}

message("SPARK-X benchmark complete.")
