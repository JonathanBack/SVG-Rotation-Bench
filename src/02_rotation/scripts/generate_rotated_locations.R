project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

input_location_file <- file.path(project_root, "src", "01_simulation", "outputs", "scDesign3", "data", "location.csv")
rotation_output_dir <- file.path(project_root, "src", "02_rotation", "outputs", "locations")
dir.create(rotation_output_dir, recursive = TRUE, showWarnings = FALSE)

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

rotate_coordinates <- function(coordinates, angle_degrees, center = TRUE) {
  coordinates <- as.matrix(coordinates)
  if (ncol(coordinates) != 2) {
    stop("`coordinates` must have exactly two columns.", call. = FALSE)
  }
  if (is.null(colnames(coordinates))) {
    colnames(coordinates) <- c("spatial1", "spatial2")
  }
  center_point <- if (center) {
    colMeans(coordinates, na.rm = TRUE)
  } else {
    c(0, 0)
  }
  translated <- sweep(coordinates, 2, center_point, FUN = "-")
  rotated <- translated %*% t(rotation_matrix_2d(angle_degrees))
  rotated <- sweep(rotated, 2, center_point, FUN = "+")
  colnames(rotated) <- colnames(coordinates)
  rownames(rotated) <- rownames(coordinates)
  rotated
}

location_table <- read.csv(input_location_file, row.names = 1, check.names = FALSE)
angles_degrees <- c(0, 30, 45, 60)

coords <- as.matrix(location_table[, c("spatial1", "spatial2")])
metadata <- location_table[, setdiff(names(location_table), c("spatial1", "spatial2")), drop = FALSE]

for (angle in angles_degrees) {
  rot_coords <- rotate_coordinates(coords, angle_degrees = angle, center = TRUE)
  rot_df <- data.frame(
    spot_id = rownames(location_table),
    metadata,
    spatial1 = unname(rot_coords[, "spatial1"]),
    spatial2 = unname(rot_coords[, "spatial2"]),
    row.names = NULL,
    check.names = FALSE
  )
  write.csv(rot_df, file.path(rotation_output_dir, paste0("rotated_locations_", angle, ".csv")), row.names = FALSE)
}
