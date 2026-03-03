# dist_perpendicular -----------------------------------------------------------

test_that("dist_perpendicular returns 0 for identical segments", {
  expect_equal(dist_perpendicular(c(0,0), c(10,0), c(0,0), c(10,0)), 0)
})

test_that("dist_perpendicular returns 0 for collinear segments", {
  # no lateral offset at all -> both perp distances are 0
  expect_equal(dist_perpendicular(c(0,0), c(5,0), c(7,0), c(10,0)), 0)
})

test_that("dist_perpendicular is correct for uniform parallel offset", {
  # Lj sits 3 units above Li, both horizontal -> both projections land on Li
  expect_equal(dist_perpendicular(c(0,0), c(10,0), c(2,3), c(8,3)), 3)
})

test_that("dist_perpendicular uses Lehmer mean, not arithmetic mean", {
  # l_perp1 = 1, l_perp2 = 3
  # arithmetic mean = 2, Lehmer mean = (1+9)/(1+3) = 2.5
  si <- c(0,0); ei <- c(10,0)
  sj <- c(2,1); ej <- c(8,3)
  expect_equal(dist_perpendicular(si, ei, sj, ej), 2.5)
  expect_gt(dist_perpendicular(si, ei, sj, ej), 2)
})

test_that("dist_perpendicular projects onto infinite line, not just Li segment", {
  # Lj lies beyond ei on x-axis, offset 2 units vertically
  # projection falls outside Li but perp distance is still well-defined
  si <- c(0,0); ei <- c(5,0)
  sj <- c(7,2); ej <- c(9,2)
  expect_equal(dist_perpendicular(si, ei, sj, ej), 2)
})

test_that("dist_perpendicular handles degenerate Li (point) without error", {
  # fallback: Euclidean dist from si to sj = sqrt(3^2 + 4^2) = 5
  expect_equal(dist_perpendicular(c(0,0), c(0,0), c(3,4), c(8,8)), 5)
})

test_that("dist_perpendicular works for non-horizontal Li", {
  # Li von (0,0) nach (5,5), Lj parallel mit Offset
  expect_equal(dist_perpendicular(c(0,0), c(5,5), c(1,0), c(6,5)),
               sqrt(0.5), tolerance = 1e-10)
})

# dist_parallel ----------------------------------------------------------------

test_that("dist_parallel is correct when Lj lies inside Li", {
  # ps=(2,0): min(2,8)=2 | pe=(8,0): min(8,2)=2 -> min(2,2) = 2
  expect_equal(dist_parallel(c(0,0), c(10,0), c(2,0), c(8,0)), 2)
})

test_that("dist_parallel when Lj overshoots Li on the left", {
  # ps=(-3,0): min(3,13)=3 | pe=(8,0): min(8,2)=2 -> min(3,2) = 2
  expect_equal(dist_parallel(c(0,0), c(10,0), c(-3,0), c(8,0)), 2)
})

test_that("dist_parallel when Lj lies completely outside Li", {
  # both projections fall past ei=(10,0) at x=12 and x=15
  # l_par1 = min(12,2)=2, l_par2 = min(15,5)=5 -> min = 2
  expect_equal(dist_parallel(c(0,0), c(10,0), c(12,0), c(15,0)), 2)
})

test_that("dist_parallel is 0 when Lj shares endpoint with Li", {
  # pe lands exactly on ei -> l_par2 = 0 -> result = 0
  expect_equal(dist_parallel(c(0,0), c(10,0), c(5,0), c(10,0)), 0)
})

test_that("dist_parallel handles degenerate Li (point) without error", {
  expect_equal(dist_parallel(c(0,0), c(0,0), c(3,4), c(8,8)), 0)
})


# dist_angle -------------------------------------------------------------------

test_that("dist_angle returns 0 for identical direction", {
  expect_equal(dist_angle(c(0,0), c(10,0), c(0,0), c(10,0)), 0)
})

test_that("dist_angle returns lj_len for perpendicular segments", {
  # theta = 90 degrees -> full length of Lj as penalty
  expect_equal(dist_angle(c(0,0), c(10,0), c(0,0), c(0,5)), 5)
})

test_that("dist_angle scales correctly for 30 degree offset", {
  # Lj at 30 deg to x-axis, length 6 -> 6 * sin(pi/6) = 3
  lj_len <- 6
  ej <- lj_len * c(cos(pi/6), sin(pi/6))
  expect_equal(dist_angle(c(0,0), c(10,0), c(0,0), ej),
               lj_len * sin(pi/6), tolerance = 1e-10)
})

test_that("dist_angle returns lj_len for antiparallel segments", {
  # theta = 180 deg >= 90 -> full penalty
  expect_equal(dist_angle(c(0,0), c(10,0), c(10,0), c(0,0)), 10)
})

test_that("dist_angle returns 0 for degenerate Lj (point)", {
  expect_equal(dist_angle(c(0,0), c(10,0), c(3,3), c(3,3)), 0)
})


# dist_segments ----------------------------------------------------------------

test_that("dist_segments returns 0 for identical segments", {
  expect_equal(dist_segments(c(0,0), c(10,0), c(0,0), c(10,0)), 0)
})

test_that("dist_segments is symmetric", {
  si <- c(0,0); ei <- c(10,0)
  sj <- c(2,3); ej <- c(9,3)
  expect_equal(dist_segments(si, ei, sj, ej),
               dist_segments(sj, ej, si, ei))
})

test_that("dist_segments swaps Li/Lj so the longer segment is always Li", {
  # Call with short Li and long Lj; result must match the correctly ordered call
  si_short <- c(0,0); ei_short <- c(4,0)
  si_long  <- c(0,0); ei_long  <- c(10,0)
  expect_equal(dist_segments(si_short, ei_short, si_long, ei_long),
               dist_segments(si_long,  ei_long,  si_short, ei_short))
})

test_that("dist_segments respects custom weights", {
  si <- c(0,0); ei <- c(10,0)
  # Lj is tilted (not parallel to Li) so dist_angle > 0
  sj <- c(0,0); ej <- c(8,3)
  d_full     <- dist_segments(si, ei, sj, ej)
  d_no_angle <- dist_segments(si, ei, sj, ej, w_angle = 0)
  d_no_perp  <- dist_segments(si, ei, sj, ej, w_perp  = 0)
  expect_gt(d_full, d_no_angle)
  expect_gt(d_full, d_no_perp)
})

test_that("dist_segments returns 0 when all weights are zero", {
  expect_equal(
    dist_segments(c(0,0), c(10,0), c(2,3), c(8,3),
                  w_perp = 0, w_par = 0, w_angle = 0),
    0
  )
})
