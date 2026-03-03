#' Build sparse neighbourhood list via Rcpp
#'
#' Calls \code{compute_neighbourhoods}, an Rcpp-compiled C++ function defined in \code{src/traclus_rcpp.cpp}.
#' One C++ call computes all pairwise distances and returns only the neighbour indices.
#'
#' @param segs Data frame from \code{partition_trajectories}.
#' @param eps,w_perp,w_par,w_angle Passed to \code{compute_neighbourhoods}.
#' @return A list of length n; each element is an integer vector of
#'   1-based neighbour indices.
#' @keywords internal
.build_neighbourhood_list <- function(segs, eps, w_perp, w_par, w_angle) {
  compute_neighbourhoods(
    sx = segs$sx, sy = segs$sy,
    ex = segs$ex, ey = segs$ey,
    eps     = eps,
    w_perp  = w_perp,
    w_par   = w_par,
    w_angle = w_angle
  )
}

#' Extract neighbours from precomputed list
#'
#' @param idx Integer row index of the query segment (1-based).
#' @param nbr_list List from \code{.build_neighbourhood_list}.
#' @return Integer vector of neighbour indices.
#' @keywords internal
.eps_neighbourhood <- function(idx, nbr_list) {
  nbr_list[[idx]]
}


#' Density-based clustering of line segments (TRACLUS grouping phase)
#'
#' Adapts the DBSCAN algorithm to line segments using the TRACLUS distance function.
#' A segment qualifies as a core segment when its epsilon-neighbourhood contains at least
#' \code{min_lns} segments. After clustering, a trajectory cardinality check discards any cluster
#' whose segments originate from fewer than \code{min_lns} distinct trajectories.
#'
#' @param segs Data frame produced by \code{partition_trajectories}, with
#'   columns \code{traj_id}, \code{sx}, \code{sy}, \code{ex}, \code{ey}.
#' @param eps Non-negative distance threshold for neighbourhood queries.
#' @param min_lns Minimum number of segments in an epsilon-neighbourhood for a segment to qualify as a
#' core segment. Also used as the minimum trajectory cardinality threshold for cluster retention.
#' @param w_perp,w_par,w_angle Non-negative weights passed to \code{dist_segments}. All default to 1.
#' @return A copy of \code{segs} with an additional integer column
#'   \code{cluster_id}. Noise segments are assigned \code{NA}.
#' @export
#' @examples
#' trajs  <- make_toy_trajectories()
#' segs   <- partition_trajectories(trajs)
#' result <- cluster_segments(segs, eps = 10, min_lns = 3)
#' table(result$cluster_id, useNA = "ifany")
cluster_segments <- function(segs,
                             eps,
                             min_lns,
                             w_perp  = 1,
                             w_par   = 1,
                             w_angle = 1) {
  stopifnot(
    is.data.frame(segs),
    all(c("traj_id", "sx", "sy", "ex", "ey") %in% names(segs)),
    eps     >= 0,
    min_lns >= 1L
  )

  n <- nrow(segs)
  # Compute all pairwise neighbourhoods in one C++ call — avoids n^2 R dispatch overhead
  nbr_mat <- .build_neighbourhood_list(segs, eps, w_perp, w_par, w_angle)

  # Cluster assignments: NA = unvisited/noise, integer = cluster ID
  cluster_id <- rep(NA_integer_, n)
  visited    <- rep(FALSE, n)
  next_id    <- 1L

  # Pre-allocate queue buffer — reused across clusters to avoid repeated c() copies
  queue <- integer(n)

  for (i in seq_len(n)) {
    # Skip segments already assigned during a previous core segment's expansion
    if (visited[i]) next
    visited[i] <- TRUE

    # Retrieve precomputed neighbours of segment i
    nbrs <- .eps_neighbourhood(i, nbr_mat)

    # Core segment check: |N_eps(L_i)| >= MinLns  (Definition 5 in Lee et al. 2007)
    if (length(nbrs) < min_lns) next

    # i is a core segment -> assign only unassigned neighbours to the new cluster
    unassigned        <- nbrs[is.na(cluster_id[nbrs])]
    cluster_id[unassigned] <- next_id
    cluster_id[i]     <- next_id

    # Seed the queue with unvisited neighbours (excluding i itself)
    seed    <- unassigned[unassigned != i]
    q_tail  <- length(seed)
    queue[seq_len(q_tail)] <- seed

    q_head <- 1L
    while (q_head <= q_tail) {
      q      <- queue[q_head]
      q_head <- q_head + 1L

      if (!visited[q]) {
        visited[q] <- TRUE
        new_nbrs   <- .eps_neighbourhood(q, nbr_mat)

        # q is also a core segment -> absorb its unassigned neighbours
        if (length(new_nbrs) >= min_lns) {
          absorb <- new_nbrs[is.na(cluster_id[new_nbrs])]
          cluster_id[absorb] <- next_id
          # Enqueue only unvisited segments to avoid redundant work
          to_enqueue <- absorb[!visited[absorb]]
          if (length(to_enqueue) > 0L) {
            queue[(q_tail + 1L):(q_tail + length(to_enqueue))] <- to_enqueue
            q_tail <- q_tail + length(to_enqueue)
          }
        }
      }

      # q may be a border segment reached from multiple cores -> assign to current cluster
      if (is.na(cluster_id[q])) cluster_id[q] <- next_id
    }

    next_id <- next_id + 1L
  }

  # Step 3 (Figure 12, lines 13-16): discard clusters whose segments originate
  # from fewer than min_lns distinct trajectories
  ids_no_na <- cluster_id[!is.na(cluster_id)]
  if (length(ids_no_na) > 0L) {
    tids_no_na <- segs$traj_id[!is.na(cluster_id)]
    n_sources  <- tapply(tids_no_na, ids_no_na, function(x) length(unique(x)))
    valid_ids  <- as.integer(names(n_sources[n_sources >= min_lns]))
    cluster_id[!is.na(cluster_id) & !(cluster_id %in% valid_ids)] <- NA_integer_
  }

  # Close gaps in cluster IDs left by discarded clusters (e.g. 1, 3, 5 -> 1, 2, 3)
  segs$cluster_id <- .renumber_clusters(cluster_id)
  segs
}


#' Re-number cluster IDs to a compact 1..K sequence
#'
#' After discarding low-cardinality clusters, the remaining IDs may be non-consecutive. This resets
#' them to 1, 2, ..., K while preserving NA for noise.
#'
#' @param ids Integer vector of cluster assignments (NA = noise).
#' @return Integer vector of the same length with reassigned IDs.
#' @keywords internal
.renumber_clusters <- function(ids) {
  # match() maps each id to its position in unique_ids and returns NA for NA inputs
  unique_ids <- sort(unique(ids[!is.na(ids)]))
  match(ids, unique_ids)
}
