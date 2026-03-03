#' @importFrom Rcpp sourceCpp
#' @useDynLib TRACLUS, .registration = TRUE
NULL

#' Parse the HURDAT2 file into a list of Cartesian trajectories
#'
#' Reads the NHC HURDAT2 best-track format and projects geographic coordinates
#' to a Cartesian plane in kilometers for use with TRACLUS.
#'
#' @details
#' HURDAT2 interleaves header rows (storm ID, name, observation count) with
#' data rows (fix observations) in a flat file. Headers are identified by the
#' Atlantic/Eastern Pacific/Central Pacific basin prefixes AL, EP, CP.
#'
#' Lat/lon are converted via equirectangular projection centred on the dataset
#' mean. Over the Atlantic basin (roughly 10N-60N, 20W-100W) the distance
#' error relative to the Haversine formula stays below 2%, which is acceptable
#' for clustering purposes.
#'
#' @param filepath Path to the HURDAT2 \code{.txt} file.
#' @param min_points Storms with fewer observations than this are dropped.
#'   They cannot produce meaningful line segments. Default: 3.
#'
#' @return A named list of trajectory matrices with columns \code{x}, \code{y}
#'   (km, Cartesian). Names are HURDAT2 storm IDs (e.g. \code{"AL011851"}).
#' @export
#' @examples
#' \dontrun{
#' path  <- system.file("extdata", "hurdat2-1851-2024-040425.txt",
#'                      package = "TRACLUS")
#' trajs <- read_hurdat2(path)
#' plot_trajectories(trajs[1:10])
#' }
read_hurdat2 <- function(filepath, min_points = 3L) {
  if (!file.exists(filepath))
    stop("File not found: ", filepath)

  raw <- readLines(filepath, warn = FALSE)

  result      <- list()
  current_id  <- NULL
  current_pts <- list()

  for (line in raw) {
    if (grepl("^(AL|EP|CP)", line)) {
      # Commit the previous storm before starting a new one
      if (!is.null(current_id) && length(current_pts) >= min_points) {
        result[[current_id]] <- do.call(rbind, current_pts)
      }
      # Extract the storm ID from the first field of the header row
      current_id  <- trimws(strsplit(line, ",")[[1L]][1L])
      current_pts <- list()

    } else {
      fields <- strsplit(line, ",")[[1L]]
      if (length(fields) < 6L) next

      # Fields 5 and 6 contain latitude and longitude with hemisphere suffixes
      lat_raw <- trimws(fields[5L])
      lon_raw <- trimws(fields[6L])

      lat <- as.numeric(sub("[NS]$", "", lat_raw))
      lon <- as.numeric(sub("[EW]$", "", lon_raw))

      # Southern and western hemispheres require negation
      if (grepl("S$", lat_raw)) lat <- -lat
      if (grepl("W$", lon_raw)) lon <- -lon

      if (!is.na(lat) && !is.na(lon))
        current_pts[[length(current_pts) + 1L]] <- c(lat = lat, lon = lon)
    }
  }

  # The file ends without a trailing header so the last storm must be flushed manually
  if (!is.null(current_id) && length(current_pts) >= min_points)
    result[[current_id]] <- do.call(rbind, current_pts)

  result <- .project_to_cartesian(result)
  message(sprintf("Loaded %d trajectories", length(result)))
  result
}


#' Project lat/lon trajectory list to Cartesian coordinates in km
#'
#' Equirectangular projection with the dataset centroid as the reference point.
#' Using raw degree differences as coordinates would make x and y distances
#' incomparable at mid-latitudes: 1 degree longitude at 50N is ~71 km,
#' while 1 degree latitude is always ~111 km.
#'
#' @param trajectories Named list of matrices with columns \code{lat}, \code{lon}
#'   (decimal degrees).
#' @return Named list of matrices with columns \code{x}, \code{y} (km).
#' @keywords internal
.project_to_cartesian <- function(trajectories) {
  R_EARTH <- 6371.0

  # Collect all coordinates to compute the dataset centroid as the projection origin
  all_lat <- unlist(lapply(trajectories, function(t) t[, "lat"]))
  all_lon <- unlist(lapply(trajectories, function(t) t[, "lon"]))

  # Reference point in radians: origin of the local Cartesian frame
  lat_ref <- mean(all_lat) * pi / 180
  lon_ref <- mean(all_lon) * pi / 180

  # Project each trajectory to the local Cartesian frame
  projected <- lapply(trajectories, function(traj) {
    lat_rad <- traj[, "lat"] * pi / 180
    lon_rad <- traj[, "lon"] * pi / 180
    # x scales longitude differences by cos(lat_ref) to correct for meridian convergence
    x <- R_EARTH * (lon_rad - lon_ref) * cos(lat_ref)
    # y is a direct arc length along the meridian
    y <- R_EARTH * (lat_rad - lat_ref)
    matrix(c(x, y), ncol = 2L,
           dimnames = list(NULL, c("x", "y")))
  })

  # Attach the projection origin so downstream code can back-project to lat/lon
  attr(projected, "lat_ref") <- lat_ref
  attr(projected, "lon_ref") <- lon_ref
  projected
}
