# Helper: straight line trajectory — should produce exactly one segment
straight_traj <- function() {
  matrix(c(0,0, 5,0, 10,0, 15,0), ncol = 2, byrow = TRUE,
         dimnames = list(NULL, c("x", "y")))
}

# Helper: L-shaped trajectory — the corner must become a characteristic point
bent_traj <- function() {
  matrix(c(0,0, 5,0, 10,0, 10,5, 10,10), ncol = 2, byrow = TRUE,
         dimnames = list(NULL, c("x", "y")))
}

test_that("straight trajectory produces one segment", {
  segs <- partition_trajectory(straight_traj())
  expect_equal(length(segs), 1L)
})

test_that("segment endpoints match trajectory start and end", {
  traj <- straight_traj()
  segs <- partition_trajectory(traj)
  expect_equal(segs[[1]]$si, traj[1,  ])
  expect_equal(segs[[length(segs)]]$ei, traj[nrow(traj), ])
})

test_that("bent trajectory produces at least two segments", {
  segs <- partition_trajectory(bent_traj())
  expect_gte(length(segs), 2L)
})

test_that("consecutive segments share endpoints", {
  segs <- partition_trajectory(bent_traj())
  for (k in seq_len(length(segs) - 1L)) {
    expect_equal(segs[[k]]$ei, segs[[k + 1L]]$si)
  }
})

test_that("two-point trajectory produces exactly one segment", {
  traj <- matrix(c(0,0, 10,5), ncol = 2, byrow = TRUE,
                 dimnames = list(NULL, c("x", "y")))
  segs <- partition_trajectory(traj)
  expect_equal(length(segs), 1L)
})

test_that("partition_trajectories returns a data frame with correct columns", {
  result <- partition_trajectories(make_toy_trajectories())
  expect_s3_class(result, "data.frame")
  expect_named(result, c("traj_id", "seg_id", "sx", "sy", "ex", "ey"))
})

test_that("partition_trajectories preserves trajectory names", {
  result <- partition_trajectories(make_toy_trajectories())
  expect_true(all(c("TR1", "TR2", "TR3", "TR4", "TR5") %in% result$traj_id))
})
