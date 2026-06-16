project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

source(file.path(project_root, "src", "02_rotation", "rotation_utils.R"))

input_location_file <- file.path(project_root, "src", "01_simulation", "outputs", "scDesign3", "location.csv")
rotation_output_dir <- file.path(project_root, "src", "02_rotation", "outputs")
rotation_file <- file.path(rotation_output_dir, "rotated_locations.csv")

dir.create(rotation_output_dir, recursive = TRUE, showWarnings = FALSE)

location_table <- read.csv(input_location_file, row.names = 1, check.names = FALSE)
angles_degrees <- c(0, 30, 45, 60)

rotated_locations <- build_rotated_location_table(
  location_table = location_table,
  angles_degrees = angles_degrees,
  center = TRUE
)

write.csv(rotated_locations, rotation_file, row.names = FALSE)

rotation_manifest <- lapply(split(rotated_locations, rotated_locations$angle_degrees), function(angle_table) {
  angle_value <- unique(angle_table$angle_degrees)
  angle_file <- file.path(rotation_output_dir, paste0("rotated_locations_", angle_value, ".csv"))

  write.csv(angle_table, angle_file, row.names = FALSE)

  data.frame(
    angle_degrees = angle_value,
    file = angle_file,
    n_spots = nrow(angle_table),
    stringsAsFactors = FALSE
  )
})

write.csv(do.call(rbind, rotation_manifest), file.path(rotation_output_dir, "rotation_manifest.csv"), row.names = FALSE)

invisible(rotated_locations)
