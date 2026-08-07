# ==============================================================================
# generate_rotated_locations.R
# Applies 2D rotation (0, 30, 45, 60 degrees) to the spatial coordinates of
# each (slice, mode) dataset. Each rotated set is saved as a separate CSV.
#
# Edit SLICES and MODES below to control scope.
# Input:  src/01_simulation/outputs/scDesign3/{slice}/{mode}/data/location.csv
# Output: src/02_rotation/outputs/{slice}/{mode}/locations/rotated_locations_{angle}.csv
# ==============================================================================

# --- Select which slices and modes to run ---
SLICES <- c("anterior1", "anterior2", "posterior1", "posterior2")
MODES  <- c("simulated", "whole")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)



angles_degrees <- c(0, 30, 45, 60)

for (slice in SLICES) {
  for (mode in MODES) {

    input_location_file <- file.path(project_root, "src", "01_simulation",
                                     "outputs", "scDesign3", slice, mode, "data",
                                     "location.csv")
    rotation_output_dir <- file.path(project_root, "src", "02_rotation",
                                     "outputs", slice, mode, "locations")
    dir.create(rotation_output_dir, recursive = TRUE, showWarnings = FALSE)

    if (!file.exists(input_location_file)) {
      message("Skipping ", slice, "/", mode, " -- location.csv not found")
      next
    }
    message("Rotating coordinates for: ", slice, " / ", mode)

    # --- Load original locations and extract coordinate columns ---
    location_table <- read.csv(input_location_file, row.names = 1, check.names = FALSE)

    coords <- as.matrix(location_table[, c("spatial1", "spatial2")])
    # Preserve non-coordinate metadata columns (e.g., cell_type)
    metadata <- location_table[, setdiff(names(location_table),
                                          c("spatial1", "spatial2")),
                                drop = FALSE]
    if (ncol(metadata) == 0) metadata <- NULL

    # --- Rotate and export for each angle ---
    for (angle in angles_degrees) {
      angle_radians <- angle * pi / 180
      R <- matrix(
        c(cos(angle_radians),  sin(angle_radians),
          -sin(angle_radians), cos(angle_radians)),
        nrow = 2, byrow = TRUE
      )
      rot_coords <- as.matrix(coords) %*% R
      colnames(rot_coords) <- c("spatial1", "spatial2")
      rownames(rot_coords) <- rownames(coords)

      rot_df <- data.frame(
        spot_id = rownames(location_table),
        spatial1 = unname(rot_coords[, "spatial1"]),
        spatial2 = unname(rot_coords[, "spatial2"]),
        row.names = NULL,
        check.names = FALSE
      )
      if (!is.null(metadata) && ncol(metadata) > 0) {
        rot_df <- cbind(rot_df, metadata)
      }
      write.csv(rot_df,
                file.path(rotation_output_dir,
                           paste0("rotated_locations_", angle, ".csv")),
                row.names = FALSE)
    }

    message("  Done: ", slice, " / ", mode)
  }
}

message("\nAll rotated location files generated.")