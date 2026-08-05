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

# --- Rotation matrix builder: standard 2D rotation by angle_degrees ---
rotation_matrix_2d <- function(angle_degrees) {
  angle_radians <- angle_degrees * pi / 180
  matrix(
    c(
      cos(angle_radians), -sin(angle_radians),
      sin(angle_radians),  cos(angle_radians)
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("spatial1", "spatial2"), c("spatial1", "spatial2"))
  )
}

# --- Coordinate rotation: center, rotate, uncenter ---
#   coordinates: n x 2 matrix (spatial1, spatial2)
#   angle_degrees: rotation angle in degrees
#   center: if TRUE, rotation is performed around the centroid (default)
rotate_coordinates <- function(coordinates, angle_degrees, center = TRUE) {
  coordinates <- as.matrix(coordinates)
  if (ncol(coordinates) != 2) {
    stop("`coordinates` must have exactly two columns.", call. = FALSE)
  }
  if (is.null(colnames(coordinates))) {
    colnames(coordinates) <- c("spatial1", "spatial2")
  }
  # Compute centroid (if centering) or use origin
  center_point <- if (center) {
    colMeans(coordinates, na.rm = TRUE)
  } else {
    c(0, 0)
  }
  # Translate to origin, apply rotation, translate back
  translated <- sweep(coordinates, 2, center_point, FUN = "-")
  rotated <- translated %*% t(rotation_matrix_2d(angle_degrees))
  rotated <- sweep(rotated, 2, center_point, FUN = "+")
  colnames(rotated) <- colnames(coordinates)
  rownames(rotated) <- rownames(coordinates)
  rotated
}

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
      rot_coords <- rotate_coordinates(coords, angle_degrees = angle, center = FALSE)
      rot_df <- data.frame(
        spot_id = rownames(location_table),
        metadata,
        spatial1 = unname(rot_coords[, "spatial1"]),
        spatial2 = unname(rot_coords[, "spatial2"]),
        row.names = NULL,
        check.names = FALSE
      )
      write.csv(rot_df,
                file.path(rotation_output_dir,
                           paste0("rotated_locations_", angle, ".csv")),
                row.names = FALSE)
    }

    message("  Done: ", slice, " / ", mode)
  }
}

message("\nAll rotated location files generated.")