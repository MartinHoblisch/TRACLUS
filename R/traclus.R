#' @importFrom stats quantile
NULL

#' Trajectory clustering via the partition-and-group framework
#'
#' Full TRACLUS pipeline: partitions each trajectory into line segments using the MDL principle, clusters
#' the segments with a density-based approach, and computes a representative trajectory for each cluster.
#'
#' @param trajectories A named list of trajectory matrices. Each matrix must have columns \code{x}
#' and \code{y} (Cartesian coordinates, e.g. km).
#'   Use \code{read_hurdat2} to obtain a correctly projected list from raw HURDAT2 data.
#' @param eps Non-negative distance threshold for the neighbourhood query.
#'   Segments within \code{eps} of each other are considered neighbours.
#'   Use \code{estimate_eps} to get a data-driven starting value.
#' @param min_lns Minimum neighbourhood size for a segment to be a core segment. Also controls the
#' minimum number of distinct source trajectories a cluster must span to be retained.
#' @param w_perp,w_par,w_angle Non-negative weights for the three distance components.
#' @param gamma Minimum spacing (same unit as coordinates) between consecutive waypoints of a
#' representative trajectory. Default 1.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{segments}}{Data frame of all line segments with a
#'       \code{cluster_id} column (\code{NA} = noise).}
#'     \item{\code{representatives}}{Named list of representative trajectory matrices (one per cluster).}
#'     \item{\code{params}}{List of the parameters used for this run.}
#'   }
#' @export
#' @examples
#' trajs  <- make_toy_trajectories()
#' result <- traclus(trajs, eps = 10, min_lns = 3)
#' length(result$representatives)
#' table(result$segments$cluster_id, useNA = "ifany")
traclus <- function(trajectories,
                    eps,
                    min_lns,
                    w_perp  = 1,
                    w_par   = 1,
                    w_angle = 1,
                    gamma   = 1) {
  stopifnot(
    is.list(trajectories), length(trajectories) >= 1L,
    eps     >= 0,
    min_lns >= 1L,
    w_perp  >= 0, w_par >= 0, w_angle >= 0,
    gamma   >= 0
  )

  message(sprintf("Partitioning %d trajectories...", length(trajectories)))
  # Phase 1: convert each trajectory into a set of characteristic line segments
  segs <- partition_trajectories(trajectories)
  message(sprintf("  -> %d line segments", nrow(segs)))

  message("Clustering segments...")
  # Phase 2: group segments by spatial and directional similarity
  clustered <- cluster_segments(segs,
                                eps     = eps,
                                min_lns = min_lns,
                                w_perp  = w_perp,
                                w_par   = w_par,
                                w_angle = w_angle)

  n_clusters <- length(unique(clustered$cluster_id[!is.na(clustered$cluster_id)]))
  n_noise    <- sum(is.na(clustered$cluster_id))
  message(sprintf("  -> %d clusters, %d noise segments", n_clusters, n_noise))

  message("Computing representative trajectories...")
  # Phase 3: summarise each cluster into a single representative path
  reps <- compute_all_representatives(clustered,
                                      min_lns = min_lns,
                                      gamma   = gamma)
  message(sprintf("  -> %d representative trajectories", length(reps)))

  # Return all three outputs together with the parameters used for reproducibility
  list(
    segments        = clustered,
    representatives = reps,
    params          = list(
      eps     = eps,
      min_lns = min_lns,
      w_perp  = w_perp,
      w_par   = w_par,
      w_angle = w_angle,
      gamma   = gamma
    )
  )
}


#' Estimate a starting value for the eps parameter
#'
#' Computes the neighbourhood-size entropy over a grid of candidate epsilon values and returns the value
#' at the entropy minimum. A skewed neighbourhood size distribution (low entropy) indicates that the
#' distance threshold separates dense regions from sparse ones effectively.
#'
#' This is a direct implementation of the heuristic in Section 4.4 of Lee et al. (2007).
#' The result is a reasonable starting point; the user should inspect a few values around it.
#'
#' @param segs Data frame produced by \code{partition_trajectories}.
#' @param eps_grid Numeric vector of candidate epsilon values to evaluate.
#'   If \code{NULL} (default), a grid of 30 values between the 5th and 95th
#'   percentile of all pairwise distances in a random subsample is used.
#' @param sample_size Number of segments to subsample for the pairwise distance
#'   estimate. Ignored if \code{eps_grid} is supplied. Default 200.
#' @param w_perp,w_par,w_angle Distance weights, passed to \code{dist_segments}.
#' @return A list with elements \code{eps_opt} (optimal epsilon) and
#'   \code{entropy_df} (data frame of eps vs entropy for diagnostic plots).
#' @export
#' @examples
#' trajs <- make_toy_trajectories()
#' segs  <- partition_trajectories(trajs)
#' est   <- estimate_eps(segs)
#' est$eps_opt
estimate_eps <- function(segs,
                         eps_grid    = NULL,
                         sample_size = 200L,
                         w_perp      = 1,
                         w_par       = 1,
                         w_angle     = 1) {
  n <- nrow(segs)

  # Always subsample: used both for automatic eps_grid construction and entropy evaluation.
  # Keeping the subsample consistent avoids comparing distances on different segment sets.
  idx <- sample(seq_len(n), min(n, sample_size))
  sub <- segs[idx, ]

  if (is.null(eps_grid)) {
    m <- nrow(sub)

    # Compute all pairwise distances on the subsample to determine a sensible eps range
    dists <- numeric(m * (m - 1L) / 2L)
    k     <- 0L
    for (i in seq_len(m - 1L)) {
      for (j in (i + 1L):m) {
        k        <- k + 1L
        dists[k] <- dist_segments(
          c(sub$sx[i], sub$sy[i]), c(sub$ex[i], sub$ey[i]),
          c(sub$sx[j], sub$sy[j]), c(sub$ex[j], sub$ey[j]),
          w_perp = w_perp, w_par = w_par, w_angle = w_angle
        )
      }
    }
    # Exclude zero distances (identical segments) and focus on the central range
    dists    <- dists[dists > 0]
    eps_grid <- seq(quantile(dists, 0.05), quantile(dists, 0.95),
                    length.out = 30L)
  }

  # Evaluate entropy of the neighbourhood size distribution at each candidate eps.
  # Low entropy means neighbourhood sizes are skewed: a few dense, most sparse.
  entropy_vals <- numeric(length(eps_grid))

  for (g in seq_along(eps_grid)) {
    eps_cand <- eps_grid[g]
    # One C++ call computes all neighbourhoods for the subsample at this eps
    nbr_list  <- .build_neighbourhood_list(sub, eps_cand, w_perp, w_par, w_angle)
    nbr_sizes <- vapply(nbr_list, length, integer(1L))

    total <- sum(nbr_sizes)
    # All segments isolated: entropy is zero by convention
    if (total == 0L) {
      entropy_vals[g] <- 0
      next
    }
    p  <- nbr_sizes / total
    p  <- p[p > 0]
    entropy_vals[g] <- -sum(p * log2(p))
  }

  # Return the eps value with the lowest entropy as the recommended starting point
  best <- which.min(entropy_vals)

  list(
    eps_opt    = eps_grid[best],
    entropy_df = data.frame(eps = eps_grid, entropy = entropy_vals)
  )
}
