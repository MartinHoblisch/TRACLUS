#' MDL cost with partitioning from pi to pj
#'
#' Computes L(H) + L(D|H) when treating the straight line between the points
#' at indices i and j in \code{pts} as the sole trajectory partition in that range.
#'
#' L(H) is the log2 length of the segment. L(D|H) sums up how much each raw sub-segment between
#' i and j deviates from the candidate segment.
#'
#' @param pts Numeric matrix with columns x, y. One row per point.
#' @param i Integer start index into \code{pts} (1-based).
#' @param j Integer end index into \code{pts} (1-based). Must satisfy i < j.
#' @return Scalar MDL cost (bits).
#' @keywords internal
.mdl_par <- function(pts, i, j) {
  # Start and end point of the candidate segment
  si <- pts[i, ]
  ei <- pts[j, ]

  # Length of the candidate segment
  seg_len <- sqrt(sum((ei - si)^2))

  # Degenerate case: segment has no length, no encoding cost
  if (seg_len < .Machine$double.eps) return(0)

  # Hypothesis cost: longer segments cost more bits to encode
  l_h <- log2(seg_len)

  # Data cost: sum deviations of each raw sub-segment from the candidate
  l_dh <- 0
  for (k in i:(j - 1L)) {
    # Start and end of the k-th raw sub-segment
    sk <- pts[k, ]
    ek <- pts[k + 1L, ]

    # Perpendicular and angle deviation from candidate segment
    d_perp  <- dist_perpendicular(si, ei, sk, ek)
    d_angle <- dist_angle(si, ei, sk, ek)

    # log2(0) is undefined, so only add non-zero deviations
    if (d_perp  > 0) l_dh <- l_dh + log2(d_perp)
    if (d_angle > 0) l_dh <- l_dh + log2(d_angle)
  }

  # Total cost = hypothesis cost + data cost
  l_h + l_dh
}


#' MDL cost without partitioning from pi to pj
#'
#' Computes L(H) when keeping every raw sub-segment between i and j.
#' L(D|H) is zero here because the raw segments perfectly encode the data.
#'
#' @inheritParams .mdl_par
#' @return Scalar MDL cost (bits).
#' @keywords internal
.mdl_nopar <- function(pts, i, j) {
  # Sum log2 lengths of each raw sub-segment
  total <- 0
  for (k in i:(j - 1L)) {
    # Length of the k-th raw sub-segment
    seg_len <- sqrt(sum((pts[k + 1L, ] - pts[k, ])^2))

    # Skip zero-length sub-segments (log2(0) undefined)
    if (seg_len > .Machine$double.eps)
      total <- total + log2(seg_len)
  }
  total
}


#' Partition a single trajectory into line segments using the MDL principle
#'
#' Greedy algorithm. Scans the trajectory point by point. When partitioning becomes cheaper than keeping
#' the raw sub-segments, the previous point is saved as a characteristic point and the scan restarts from there.
#'
#' @param traj Numeric matrix with columns x, y. At least 2 rows.
#' @return A list of line segments. Each element is a named list with si and ei (numeric vectors of length 2).
#' @export
#' @examples
#' traj <- matrix(c(0,0, 5,1, 10,0, 15,5, 20,10), ncol=2, byrow=TRUE,
#'                dimnames=list(NULL, c("x","y")))
#' segs <- partition_trajectory(traj)
#' length(segs)
partition_trajectory <- function(traj) {
  n <- nrow(traj)
  if (n < 2L) stop("A trajectory needs at least 2 points to be partitioned.")

  # Start with the first point as a characteristic point
  char_points <- 1L

  start   <- 1L  # index of the current segment start
  n_steps <- 1L  # how many steps were extended so far

  while (start + n_steps <= n) {
    # Current end point of the candidate segment
    curr <- start + n_steps

    # Compare MDL cost with and without partitioning
    cost_par   <- .mdl_par(traj, start, curr)
    cost_nopar <- .mdl_nopar(traj, start, curr)

    if (cost_par > cost_nopar) {
      # Partitioning is cheaper: save the previous point and restart from there
      char_points <- c(char_points, curr - 1L)
      start   <- curr - 1L
      n_steps <- 1L
    } else {
      # Keep extending the candidate segment
      n_steps <- n_steps + 1L
    }
  }

  # Always include the last point
  char_points <- c(char_points, n)
  char_points <- unique(char_points)

  # Build a list of segments from consecutive characteristic points
  segments <- vector("list", length(char_points) - 1L)
  for (k in seq_len(length(char_points) - 1L)) {
    segments[[k]] <- list(
      si = traj[char_points[k],      ],
      ei = traj[char_points[k + 1L], ]
    )
  }

  segments
}


#' Partition all trajectories in a list
#'
#' Calls partition_trajectory on each trajectory and returns all segments in a single data frame
#' with trajectory and segment IDs.
#'
#' @param trajectories A named list of trajectory matrices (columns x, y).
#' @return A data frame with columns:
#'   \describe{
#'     \item{traj_id}{Name of the source trajectory.}
#'     \item{seg_id}{Segment index within that trajectory.}
#'     \item{sx,sy}{Start point coordinates.}
#'     \item{ex,ey}{End point coordinates.}
#'   }
#' @export
#' @examples
#' trajs <- make_toy_trajectories()
#' segs  <- partition_trajectories(trajs)
#' nrow(segs)
partition_trajectories <- function(trajectories) {
  n_traj     <- length(trajectories)
  traj_names <- names(trajectories)

  # First pass: partition every trajectory and count the total number of segments
  all_segs   <- vector("list", n_traj)
  seg_counts <- integer(n_traj)

  for (i in seq_len(n_traj)) {
    all_segs[[i]]   <- partition_trajectory(trajectories[[i]])
    seg_counts[i]   <- length(all_segs[[i]])
  }

  total <- sum(seg_counts)

  # Pre-allocate flat vectors — avoids repeated data.frame() + rbind() overhead
  out_traj_id <- character(total)
  out_seg_id  <- integer(total)
  out_sx      <- numeric(total)
  out_sy      <- numeric(total)
  out_ex      <- numeric(total)
  out_ey      <- numeric(total)

  # Second pass: fill the pre-allocated vectors
  pos <- 1L
  for (i in seq_len(n_traj)) {
    tid  <- if (!is.null(traj_names) && !is.na(traj_names[i])) traj_names[i] else as.character(i)
    segs <- all_segs[[i]]

    for (j in seq_along(segs)) {
      out_traj_id[pos] <- tid
      out_seg_id[pos]  <- j
      out_sx[pos]      <- segs[[j]]$si[1]
      out_sy[pos]      <- segs[[j]]$si[2]
      out_ex[pos]      <- segs[[j]]$ei[1]
      out_ey[pos]      <- segs[[j]]$ei[2]
      pos              <- pos + 1L
    }
  }

  data.frame(
    traj_id = out_traj_id,
    seg_id  = out_seg_id,
    sx = out_sx, sy = out_sy,
    ex = out_ex, ey = out_ey,
    stringsAsFactors = FALSE
  )
}
