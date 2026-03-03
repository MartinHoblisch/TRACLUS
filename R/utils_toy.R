#' @importFrom graphics lines points
NULL

#' Small toy dataset for TRACLUS testing
#'
#' Six hand-crafted trajectories: five (TR1--TR5) share a common sub-trajectory through the middle section;
#' TR6 is an intentional outlier that crosses the others vertically. Coordinates are in kilometers (Cartesian).
#'
#' @return A named list of six trajectory matrices, each with columns
#'   \code{x} and \code{y} (km).
#' @export
#' @examples
#' trajs <- make_toy_trajectories()
#' plot_trajectories(trajs)
make_toy_trajectories <- function() {
  list(
    TR1 = matrix(c(
      0,  0,
      10,  2,
      20,  4,
      30,  4,
      40,  4,
      50,  8,
      60, 15
    ), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("x", "y"))),

    TR2 = matrix(c(
      5, -5,
      15,  1,
      25,  3,
      35,  5,
      45,  5,
      52, 12,
      55, 20
    ), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("x", "y"))),

    TR3 = matrix(c(
      10,  0,
      20,  3,
      30,  4,
      40,  5,
      50,  4,
      60,  0,
      70, -5
    ), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("x", "y"))),

    TR4 = matrix(c(
      8,  2,
      18,  4,
      28,  5,
      38,  6,
      48,  6,
      50, 15,
      52, 25
    ), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("x", "y"))),

    TR5 = matrix(c(
      0, 12,
      10,  8,
      20,  5,
      30,  4,
      40,  4,
      55,  2,
      70, -2
    ), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("x", "y"))),

    TR6 = matrix(c(
      35, 35,
      35, 25,
      36, 15,
      37,  5,
      38, -5
    ), ncol = 2, byrow = TRUE, dimnames = list(NULL, c("x", "y")))
  )
}


#' Quick trajectory plot for visual inspection
#'
#' Base-R only, intentionally minimal. This is a debugging aid, not a publication figure.
#'
#' @param trajectories A list of trajectory matrices with columns \code{x}, \code{y}.
#' @param title Plot title.
#' @param col Optional colour vector. If \code{NULL} (default), the Wong (2011) colorblind-safe palette is used,
#' cycling if there are more than 8 trajectories.
#' @param lwd Line width passed to \code{\link[graphics]{lines}}. Default 2.
#' @return Invisibly \code{NULL}, called for its side effect.
#' @export
#' @examples
#' plot_trajectories(make_toy_trajectories())
plot_trajectories <- function(trajectories, title = "Trajectories", col = NULL, lwd = 2) {
  all_x <- unlist(lapply(trajectories, function(t) t[, "x"]))
  all_y <- unlist(lapply(trajectories, function(t) t[, "y"]))

  plot(NULL,
       xlim = range(all_x),
       ylim = range(all_y),
       xlab = "x (km)",
       ylab = "y (km)",
       main = title,
       asp  = 1
  )

  # Wong (2011) colorblind-safe palette, cycles if more than 8 trajectories
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
           "#0072B2", "#D55E00", "#CC79A7", "#000000")

  for (i in seq_along(trajectories)) {
    line_col <- if (is.null(col)) pal[(i - 1L) %% length(pal) + 1L] else col

    traj <- trajectories[[i]]
    lines(traj[,  "x"], traj[, "y"], col = line_col, lwd = lwd)

    points(traj[, "x"], traj[, "y"], col = line_col, pch = 19, cex = 0.8 * (lwd/2))
  }

  invisible(NULL)
}
