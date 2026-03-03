# Fixture: three horizontal parallel segments — representative should be
# a horizontal line near their common y-centre
make_rep_fixture <- function() {
  data.frame(
    traj_id    = c("A", "B", "C"),
    seg_id     = c(1L, 1L, 1L),
    sx         = c(0,  0,  0),
    sy         = c(-1, 0,  1),
    ex         = c(10, 10, 10),
    ey         = c(-1, 0,  1),
    cluster_id = c(1L, 1L, 1L),
    stringsAsFactors = FALSE
  )
}

test_that("compute_representative returns a matrix with x and y columns", {
  rt <- compute_representative(make_rep_fixture(), min_lns = 3, gamma = 0)
  expect_true(is.matrix(rt))
  expect_named(as.data.frame(rt), c("x", "y"))
})

test_that("representative y is near the mean of parallel segments", {
  rt <- compute_representative(make_rep_fixture(), min_lns = 3, gamma = 0)
  expect_true(all(abs(rt[, "y"]) <= 1))
})

test_that("compute_representative returns NULL when density is too low", {
  rt <- compute_representative(make_rep_fixture(), min_lns = 5, gamma = 0)
  expect_null(rt)
})

test_that("gamma filters waypoints that are too close together", {
  rt_fine   <- compute_representative(make_rep_fixture(), min_lns = 2, gamma = 0)
  rt_coarse <- compute_representative(make_rep_fixture(), min_lns = 2, gamma = 5)
  expect_gte(nrow(rt_fine), nrow(rt_coarse))
})

test_that("compute_all_representatives returns one entry per valid cluster", {
  trajs  <- make_toy_trajectories()
  segs   <- partition_trajectories(trajs)
  result <- cluster_segments(segs, eps = 10, min_lns = 3)
  rts    <- compute_all_representatives(result, min_lns = 3)
  n_clusters <- length(unique(result$cluster_id[!is.na(result$cluster_id)]))
  expect_lte(length(rts), n_clusters)
})
