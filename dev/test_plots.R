rm(list = ls())

### hurdart import
path   <- system.file("extdata", "hurdat2-1851-2024-040425.txt", package = "TRACLUS")
hurdat <- read_hurdat2(path)

### Visualization of the original dataset
library(leaflet)
library(sf)

lat_ref <- attr(hurdat, "lat_ref")
lon_ref <- attr(hurdat, "lon_ref")

km_to_lat <- function(y) (lat_ref + y / 6371.0) * 180 / pi
km_to_lon <- function(x) (lon_ref + x / (6371.0 * cos(lat_ref))) * 180 / pi

m <- leaflet() |> addProviderTiles(providers$CartoDB.Positron)

trajs_sfc <- st_sfc(lapply(hurdat_sub, function(traj) {
  coords <- matrix(c(km_to_lon(traj[, "x"]), km_to_lat(traj[, "y"])), ncol = 2)
  st_linestring(coords)
}), crs = 4326)

trajs_sf <- st_sf(geometry = trajs_sfc)

m <- m |> addPolylines(
  data    = trajs_sf,
  color   = "#333333",
  weight  = 0.5,
  opacity = 0.3
)

m

### Subset to 1950-2004, matching the evaluation period in Lee et al. (2007)
storm_years <- as.integer(substr(names(hurdat), 5, 8))
hurdat_sub  <- hurdat[storm_years >= 1950 & storm_years <= 2004]

result <- traclus(hurdat_sub, eps = 75, min_lns = 3, w_perp = 1, w_par = 1, w_angle = 1, gamma = 200)

### Visualization of Trajectory calculation
plot_traclus_result(result$segments, result$representatives,
                    mode    = "clusters",
                    title   = "Atlantic hurricanes 1950-2004 (cluster view)",
                    seg_lwd = 0.8, noise_lwd = 0.3)


##### leaflet
{
  library(leaflet)
  library(sf)
  #> Linking to GEOS 3.13.0, GDAL 3.10.1, PROJ 9.5.1; sf_use_s2() is TRUE

  # Retrieve the projection origin (in radians) used by .project_to_cartesian
  lat_ref <- attr(hurdat, "lat_ref")
  lon_ref <- attr(hurdat, "lon_ref")

  # Back-project km coordinates to geographic lat/lon (inverse of the equirectangular projection)
  km_to_lat <- function(y) (lat_ref + y / 6371.0) * 180 / pi
  km_to_lon <- function(x) (lon_ref + x / (6371.0 * cos(lat_ref))) * 180 / pi

  # Helper: convert segment rows to a single MULTILINESTRING sf object
  segs_to_sf <- function(segs) {
    lines <- lapply(seq_len(nrow(segs)), function(i) {
      st_linestring(matrix(c(
        km_to_lon(c(segs$sx[i], segs$ex[i])),
        km_to_lat(c(segs$sy[i], segs$ey[i]))
      ), ncol = 2))
    })
    st_sfc(lines, crs = 4326)
  }

  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
           "#0072B2", "#D55E00", "#CC79A7", "#000000")
  cluster_ids <- sort(unique(hurdat_clustered$cluster_id[
    !is.na(hurdat_clustered$cluster_id)]))

  m <- leaflet() |> addProviderTiles(providers$CartoDB.Positron)

  # Noise segments in grey, dashed
  noise_h <- hurdat_clustered[is.na(hurdat_clustered$cluster_id), ]
  if (nrow(noise_h) > 0L) {
    m <- m |> addPolylines(data = segs_to_sf(noise_h),
                           color = "#999999", weight = 0.5, opacity = 0.4,
                           dashArray = "4 4")
  }

  # One layer per cluster (single addPolylines call each)
  for (k in seq_along(cluster_ids)) {
    cid   <- cluster_ids[k]
    csub  <- hurdat_clustered[!is.na(hurdat_clustered$cluster_id) &
                                hurdat_clustered$cluster_id == cid, ]
    col_k <- pal[(k - 1L) %% length(pal) + 1L]
    m <- m |> addPolylines(data = segs_to_sf(csub),
                           color = col_k, weight = 1, opacity = 0.6)
  }

  # Representative trajectories as single polylines with black outline
  for (k in seq_along(hurdat_reps)) {
    rt    <- hurdat_reps[[k]]
    col_k <- pal[(k - 1L) %% length(pal) + 1L]
    coords <- matrix(c(km_to_lon(rt[, "x"]), km_to_lat(rt[, "y"])), ncol = 2)
    line_sf <- st_sfc(st_linestring(coords), crs = 4326)
    # Black outline first, then coloured line on top
    m <- m |> addPolylines(data = line_sf,
                           color = "black", weight = 7, opacity = 1)
    m <- m |> addPolylines(data = line_sf,
                           color = col_k, weight = 5, opacity = 1,
                           popup = paste("Cluster", names(hurdat_reps)[k]))
  }

  # Legend matching cluster colours
  m <- m |> addLegend(
    position = "bottomright",
    colors   = c("#999999", pal[seq_along(hurdat_reps)]),
    labels   = c("Noise", paste("Cluster", names(hurdat_reps))),
    title    = "Clusters",
    opacity  = 1
  )

  m
}


############### Toy Dataset visualize
#trajs <- make_toy_trajectories()
#segs <- partition_trajectories(trajs)
# clustered_segs <- cluster_segments(segs,
#                                    eps     = 10,
#                                    min_lns = 2,
#                                    w_perp  = 1,
#                                    w_par   = 1,
#                                    w_angle = 1)
# reps_toy <- compute_all_representatives(clustered_segs,
#                                         min_lns = 2,
#                                         gamma   = 0.1)

traclus_result <- traclus(make_toy_trajectories(),
                          eps     = 10,
                          min_lns = 2,
                          w_perp  = 1,
                          w_par   = 1,
                          w_angle = 1,
                          gamma   = 0.5)

plot_traclus_result(clustered_segs, reps_toy, title = "Clustering result (eps=20, min_lns=3)")
