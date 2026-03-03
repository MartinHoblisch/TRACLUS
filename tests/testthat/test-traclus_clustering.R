# Two tight clusters of horizontal segments, clearly separated vertically.
# Cluster A sits near y=0, cluster B near y=20. With a small eps the
# algorithm must find exactly two clusters of trajectory-diverse segments.
make_cluster_fixture <- function() {
  list(
    # Cluster A — near y = 0
    A1 = matrix(c(0,0,  10,0),  ncol=2, byrow=TRUE, dimnames=list(NULL, c("x","y"))),
    A2 = matrix(c(0,1,  10,1),  ncol=2, byrow=TRUE, dimnames=list(NULL, c("x","y"))),
    A3 = matrix(c(0,-1, 10,-1), ncol=2, byrow=TRUE, dimnames=list(NULL, c("x","y"))),
    # Cluster B — near y = 20
    B1 = matrix(c(0,20, 10,20), ncol=2, byrow=TRUE, dimnames=list(NULL, c("x","y"))),
    B2 = matrix(c(0,21, 10,21), ncol=2, byrow=TRUE, dimnames=list(NULL, c("x","y"))),
    B3 = matrix(c(0,19, 10,19), ncol=2, byrow=TRUE, dimnames=list(NULL, c("x","y")))
  )
}

test_that("cluster_segments finds two distinct clusters in fixture", {
  segs   <- partition_trajectories(make_cluster_fixture())
  result <- cluster_segments(segs, eps = 3, min_lns = 3)
  n_clusters <- length(unique(result$cluster_id[!is.na(result$cluster_id)]))
  expect_equal(n_clusters, 2L)
})

test_that("cluster IDs are re-numbered from 1", {
  segs   <- partition_trajectories(make_cluster_fixture())
  result <- cluster_segments(segs, eps = 3, min_lns = 3)
  valid  <- sort(unique(result$cluster_id[!is.na(result$cluster_id)]))
  expect_equal(valid, seq_along(valid))
})

test_that("output contains cluster_id column", {
  segs   <- partition_trajectories(make_toy_trajectories())
  result <- cluster_segments(segs, eps = 10, min_lns = 3)
  expect_true("cluster_id" %in% names(result))
})

test_that("noise segments are NA", {
  segs   <- partition_trajectories(make_cluster_fixture())
  result <- cluster_segments(segs, eps = 3, min_lns = 3)
  # cluster_id is either a positive integer or NA — nothing else
  valid  <- result$cluster_id[!is.na(result$cluster_id)]
  expect_true(all(valid >= 1L))
})

test_that("increasing eps merges or grows clusters", {
  segs      <- partition_trajectories(make_cluster_fixture())
  result_tight <- cluster_segments(segs, eps =  3, min_lns = 3)
  result_wide  <- cluster_segments(segs, eps = 25, min_lns = 3)
  n_tight <- sum(!is.na(result_tight$cluster_id))
  n_wide  <- sum(!is.na(result_wide$cluster_id))
  expect_gte(n_wide, n_tight)
})

test_that("cluster_segments stops on malformed input", {
  expect_error(cluster_segments(data.frame(x = 1), eps = 1, min_lns = 2))
})
