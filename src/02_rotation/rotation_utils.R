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

build_rotated_location_table <- function(location_table, angles_degrees, center = TRUE) {
  location_table <- as.data.frame(location_table)

  required_columns <- c("spatial1", "spatial2")
  missing_columns <- setdiff(required_columns, names(location_table))
  if (length(missing_columns) > 0) {
    stop(
      "`location_table` is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(rownames(location_table))) {
    rownames(location_table) <- seq_len(nrow(location_table))
  }

  metadata_columns <- setdiff(names(location_table), required_columns)
  metadata <- location_table[, metadata_columns, drop = FALSE]
  coordinates <- location_table[, required_columns, drop = FALSE]

  rotated_tables <- lapply(angles_degrees, function(angle_degrees) {
    rotated_coordinates <- rotate_coordinates(coordinates, angle_degrees = angle_degrees, center = center)

    data.frame(
      spot_id = rownames(location_table),
      angle_degrees = angle_degrees,
      metadata,
      spatial1 = unname(rotated_coordinates[, "spatial1"]),
      spatial2 = unname(rotated_coordinates[, "spatial2"]),
      row.names = NULL,
      check.names = FALSE
    )
  })

  do.call(rbind, rotated_tables)
}
