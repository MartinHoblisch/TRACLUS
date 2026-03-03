test_that("traclus returns a list with the three expected elements", {
  result <- traclus(make_toy_trajectories(), eps = 10, min_lns = 3)
  expect_named(result, c("segments", "representatives", "params"))
})

test_that("traclus segments data frame has cluster_id column", {
  result <- traclus(make_toy_trajectories(), eps = 10, min_lns = 3)
  expect_true("cluster_id" %in% names(result$segments))
})

test_that("traclus params match the arguments passed", {
  result <- traclus(make_toy_trajectories(),
                    eps = 7, min_lns = 2, w_perp = 2, gamma = 0.5)
  expect_equal(result$params$eps,     7)
  expect_equal(result$params$min_lns, 2)
  expect_equal(result$params$w_perp,  2)
  expect_equal(result$params$gamma,   0.5)
})

test_that("traclus finds at least one cluster on the toy dataset", {
  result <- traclus(make_toy_trajectories(), eps = 10, min_lns = 3)
  expect_gte(length(result$representatives), 1L)
})

test_that("traclus stops on invalid input", {
  expect_error(traclus(list(), eps = 10, min_lns = 3))
  expect_error(traclus(make_toy_trajectories(), eps = -1, min_lns = 3))
})

test_that("estimate_eps returns a list with eps_opt and entropy_df", {
  segs <- partition_trajectories(make_toy_trajectories())
  est  <- estimate_eps(segs)
  expect_named(est, c("eps_opt", "entropy_df"))
  expect_true(est$eps_opt > 0)
})
