#' Compute the average direction vector of a set of line segments
#'
#' Sums the raw displacement vectors (not unit vectors), so longer segments contribute proportionally
#'  more to the result.
#'
#' @param segs Data frame with columns \code{sx}, \code{sy}, \code{ex}, \code{ey}.
#' @return Numeric vector of length 2 (normalised to unit length).
#' @keywords internal
.average_direction <- function(segs) {
  # Raw displacement vectors: longer segments naturally dominate the mean direction
  dx <- segs$ex - segs$sx
  dy <- segs$ey - segs$sy

  mean_vec <- c(mean(dx), mean(dy))
  len      <- sqrt(sum(mean_vec^2))

  # Degenerate cluster: all segments cancel out, fall back to positive x-axis
  if (len < .Machine$double.eps) return(c(1, 0))
  # Normalise to unit length so the rotation matrix is orthogonal
  mean_vec / len
}


#' Rotate 2-D points so the x-axis aligns with a given direction vector
#'
#' @param pts Two-column numeric matrix (x, y).
#' @param dir Unit direction vector of length 2.
#' @return Rotated matrix with the same dimensions.
#' @keywords internal
.rotate_to_axis <- function(pts, dir) {
  # Rotation angle is implicitly encoded in dir = (cos_phi, sin_phi)
  cos_phi <- dir[1]
  sin_phi <- dir[2]
  # Standard 2-D rotation matrix applied row-wise: x' = x*cos + y*sin, y' = -x*sin + y*cos
  cbind(
    pts[, 1] * cos_phi + pts[, 2] * sin_phi,
    -pts[, 1] * sin_phi + pts[, 2] * cos_phi
  )
}


#' Rotate 2-D points back from axis-aligned frame to original frame
#'
#' @param pts Two-column numeric matrix in the rotated frame.
#' @param dir Unit direction vector used in \code{.rotate_to_axis}.
#' @return Matrix in original coordinates.
#' @keywords internal
.rotate_from_axis <- function(pts, dir) {
  cos_phi <- dir[1]
  sin_phi <- dir[2]
  # Inverse rotation: transpose of the original rotation matrix (R^-1 = R^T for orthogonal R)
  cbind(
    pts[, 1] * cos_phi - pts[, 2] * sin_phi,
    pts[, 1] * sin_phi + pts[, 2] * cos_phi
  )
}


#' Generate the representative trajectory for one cluster
#'
#' Implements the sweep-line approach: The cluster's segments are projected onto the axis of the average
#' direction vector. A sweep line then moves along that axis, wherever at least \code{min_lns} segments
#' are simultaneously active, the average transverse coordinate is recorded as a waypoint.
#'
#' @param cluster_segs Data frame of segments belonging to a single cluster
#'   (subset of the output of \code{cluster_segments}).
#' @param min_lns Minimum number of simultaneously active segments required to emit a waypoint.
#' @param gamma Minimum spacing (in km) between consecutive waypoints. Closer candidates are skipped to
#' smooth the output trajectory.
#' @return Numeric matrix with columns \code{x}, \code{y} (km). Returns
#'   \code{NULL} if no waypoints pass the density threshold.
#' @export
#' @examples
#' trajs  <- make_toy_trajectories()
#' segs   <- partition_trajectories(trajs)
#' result <- cluster_segments(segs, eps = 10, min_lns = 3)
#' csegs  <- result[!is.na(result$cluster_id) & result$cluster_id == 1, ]
#' rt     <- compute_representative(csegs, min_lns = 3, gamma = 1)
compute_representative <- function(cluster_segs, min_lns, gamma = 1) {
  stopifnot(nrow(cluster_segs) >= 1L)

  # Average direction defines the sweep axis, longer segments weighted more heavily
  dir <- .average_direction(cluster_segs)

  # Project all segment endpoints into the axis-aligned frame
  starts_rot <- .rotate_to_axis(cbind(cluster_segs$sx, cluster_segs$sy), dir)
  ends_rot   <- .rotate_to_axis(cbind(cluster_segs$ex, cluster_segs$ey), dir)

  n_segs <- nrow(cluster_segs)

  # Normalise each segment so its left endpoint (smaller x') is always in starts_rot.
  # Segments running against the average direction would otherwise fire start/end
  # events in the wrong order, causing interpolation to extrapolate outside [0, 1]
  swap <- starts_rot[, 1] > ends_rot[, 1]
  if (any(swap)) {
    tmp                   <- starts_rot[swap, , drop = FALSE]
    starts_rot[swap, ]    <- ends_rot[swap, , drop = FALSE]
    ends_rot[swap, ]      <- tmp
  }

  # Build a sorted event list: each segment contributes one start and one end event
  events <- data.frame(
    x_prime  = c(starts_rot[, 1], ends_rot[, 1]),
    seg_idx  = c(seq_len(n_segs), seq_len(n_segs)),
    is_start = c(rep(TRUE, n_segs), rep(FALSE, n_segs))
  )
  # Sorting by x' ensures the sweep processes events in spatial order
  events <- events[order(events$x_prime), ]

  # Pre-allocate waypoint storage (upper bound: one waypoint per event)
  wp_x         <- numeric(nrow(events))
  wp_y         <- numeric(nrow(events))
  wp_count     <- 0L
  active       <- logical(n_segs)  # TRUE while the sweep line intersects segment k
  last_x_prime <- -Inf

  # Cache segment coordinates as plain vectors for fast indexing inside the loop
  s_x <- starts_rot[, 1];  s_y <- starts_rot[, 2]
  e_x <- ends_rot[, 1];    e_y <- ends_rot[, 2]

  for (r in seq_len(nrow(events))) {
    idx     <- events$seg_idx[r]
    x_prime <- events$x_prime[r]

    if (events$is_start[r]) {
      # START: erst aktivieren, dann prüfen
      active[idx] <- TRUE
    }
    # END: Segment bleibt noch aktiv für diesen Check

    n_active <- sum(active)

    if (n_active >= min_lns && (x_prime - last_x_prime) >= gamma) {
      act_idx <- which(active)
      seg_dx  <- e_x[act_idx] - s_x[act_idx]
      vertical <- abs(seg_dx) < .Machine$double.eps

      t <- numeric(length(act_idx))
      t[vertical]  <- 0.5
      t[!vertical] <- (x_prime - s_x[act_idx[!vertical]]) / seg_dx[!vertical]

      y_vals <- s_y[act_idx] + t * (e_y[act_idx] - s_y[act_idx])

      wp_count       <- wp_count + 1L
      wp_x[wp_count] <- x_prime
      wp_y[wp_count] <- mean(y_vals)
      last_x_prime   <- x_prime
    }

    if (!events$is_start[r]) {
      # END: erst jetzt deaktivieren
      active[idx] <- FALSE
    }
  }

  if (wp_count == 0L) return(NULL)

  # Rotate waypoints back to the original coordinate frame
  rot_pts          <- cbind(wp_x[seq_len(wp_count)], wp_y[seq_len(wp_count)])
  result           <- .rotate_from_axis(rot_pts, dir)
  colnames(result) <- c("x", "y")
  result
}


#' Generate representative trajectories for all clusters
#'
#' Wrapper that calls \code{compute_representative} for each cluster ID found in the \code{cluster_id}
#' column of \code{clustered_segs}.
#'
#' @param clustered_segs Data frame output of \code{cluster_segments}.
#' @param min_lns Passed to \code{compute_representative}.
#' @param gamma Passed to \code{compute_representative}.
#' @return A named list of representative trajectory matrices. Names are the cluster IDs as character strings.
#' Clusters that produce no waypoints are silently dropped.
#' @export
#' @examples
#' trajs  <- make_toy_trajectories()
#' segs   <- partition_trajectories(trajs)
#' result <- cluster_segments(segs, eps = 10, min_lns = 3)
#' rts    <- compute_all_representatives(result, min_lns = 3)
#' length(rts)
compute_all_representatives <- function(clustered_segs, min_lns, gamma = 1) {
  cluster_ids <- sort(unique(clustered_segs$cluster_id[
    !is.na(clustered_segs$cluster_id)]))

  reps <- list()
  for (cid in cluster_ids) {
    csegs <- clustered_segs[!is.na(clustered_segs$cluster_id) &
                              clustered_segs$cluster_id == cid, ]
    rt <- compute_representative(csegs, min_lns = min_lns, gamma = gamma)
    if (!is.null(rt)) reps[[as.character(cid)]] <- rt
  }

  # Demote clusters without representative to noise
  valid_cids <- as.integer(names(reps))
  failed     <- !is.na(clustered_segs$cluster_id) &
    !(clustered_segs$cluster_id %in% valid_cids)
  clustered_segs$cluster_id[failed] <- NA_integer_
  clustered_segs$cluster_id <- .renumber_clusters(clustered_segs$cluster_id)
  names(reps) <- as.character(seq_along(reps))

  list(
    segments        = clustered_segs,
    representatives = reps
  )
}
