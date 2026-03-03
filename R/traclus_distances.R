# Pure-R implementation of the three TRACLUS distance components and the combined segment distance function.
#
# Note: The neighbourhood computation in the clustering phase was moved to src/traclus_rcpp.cpp for
# performance reasons — the C++ code replicates the same distance logic internally. These R functions remain
# in use for MDL partitioning and eps estimation and are preserved here as a readable reference implementation.

#' Perpendicular distance between two line segments
#'
#' Computes the Lehmer mean (order 2) of the two perpendicular distances from the endpoints of \code{Lj} to
#' the infinite line through \code{Li}.
#'
#' The Lehmer mean weights the larger distance more heavily than a plain average. This choice penalizes
#' segments that diverge sharply at one end to ensure clustering stability.
#' \deqn{d_{\perp} = \frac{l_{\perp 1}^2 + l_{\perp 2}^2}{l_{\perp 1} + l_{\perp 2}}}
#'
#' @param si,ei Numeric vectors of length 2 (x, y): start and end of \code{Li}.
#' @param sj,ej Numeric vectors of length 2 (x, y): start and end of \code{Lj}.
#' @return Non-negative scalar distance.
#' @export
#' @examples
#' dist_perpendicular(c(0,0), c(10,0), c(2,3), c(8,3))  # parallel offset
#' dist_perpendicular(c(0,0), c(10,0), c(0,0), c(10,0)) # identical -> 0
dist_perpendicular <- function(si, ei, sj, ej) {
  # Direction vector and squared length of Li
  li_vec  <- ei - si
  li_len2 <- sum(li_vec^2)

  # Avoid division by zero for degenerate segments
  if (li_len2 < .Machine$double.eps) {
    return(sqrt(sum((sj - si)^2)))
  }

  # Project sj and ej onto Li, then measure how far each sits off the line
  u_s <- sum((sj - si) * li_vec) / li_len2
  u_e <- sum((ej - si) * li_vec) / li_len2

  ps <- si + u_s * li_vec
  pe <- si + u_e * li_vec

  l_perp1_sq <- sum((sj - ps)^2)
  l_perp2_sq <- sum((ej - pe)^2)

  # Lehmer mean: reduces to the regular mean when l_perp1 == l_perp2.
  # Check denominator against machine epsilon to prevent division by zero.
  denom <- sqrt(l_perp1_sq) + sqrt(l_perp2_sq)
  if (denom < .Machine$double.eps) return(0)

  (l_perp1_sq + l_perp2_sq) / denom
}


#' Parallel distance between two line segments
#'
#' Measures how far \code{Lj} extends beyond \code{Li} along \code{Li}'s own axis.
#'
#' \code{MIN} rather than \code{MAX} is used so that a slight longitudinal overhang on one side does not
#' dominate when the segments largely overlap.
#'
#' @inheritParams dist_perpendicular
#' @return Non-negative scalar distance.
#' @export
#' @examples
#' dist_parallel(c(0,0), c(10,0), c(2,0), c(8,0))   # Lj within Li -> 2
#' dist_parallel(c(0,0), c(10,0), c(-3,0), c(8,0))  # overshoot left -> 2
dist_parallel <- function(si, ei, sj, ej) {
  li_vec  <- ei - si
  li_len2 <- sum(li_vec^2)

  # Avoid division by zero
  if (li_len2 < .Machine$double.eps) return(0)

  u_s <- sum((sj - si) * li_vec) / li_len2
  u_e <- sum((ej - si) * li_vec) / li_len2

  ps <- si + u_s * li_vec
  pe <- si + u_e * li_vec

  # Distance from each projection to the nearer endpoint of Li
  l_par1 <- min(sqrt(sum((ps - si)^2)), sqrt(sum((ps - ei)^2)))
  l_par2 <- min(sqrt(sum((pe - si)^2)), sqrt(sum((pe - ei)^2)))

  min(l_par1, l_par2)
}


#' Angle distance between two line segments
#'
#' Scales the length of \code{Lj} by \code{sin(theta)}, where theta is the smaller intersecting angle between
#' the two segments.
#'
#' When theta >= 90 degrees the full length of \code{Lj} is returned — the segments run in opposite directions,
#' which is the worst-case directional mismatch.
#'
#' @inheritParams dist_perpendicular
#' @return Non-negative scalar distance.
#' @export
#' @examples
#' dist_angle(c(0,0), c(10,0), c(0,0), c(10,0))  # same direction -> 0
#' dist_angle(c(0,0), c(10,0), c(0,0), c(0,5))   # 90 degrees -> length of Lj
dist_angle <- function(si, ei, sj, ej) {
  li_vec <- ei - si
  lj_vec <- ej - sj

  li_len <- sqrt(sum(li_vec^2))
  lj_len <- sqrt(sum(lj_vec^2))

  # A zero-length segment carries no directional information
  if (li_len < .Machine$double.eps || lj_len < .Machine$double.eps) return(0)

  cos_theta <- sum(li_vec * lj_vec) / (li_len * lj_len)
  # Floating-point arithmetic can push this just outside [-1, 1]
  cos_theta <- max(-1, min(1, cos_theta))
  theta     <- acos(cos_theta)

  if (theta < pi / 2) {
    lj_len * sin(theta)
  } else {
    lj_len
  }
}


#' Combined TRACLUS distance between two line segments
#'
#' Weighted sum of the three component distances. The longer segment is always assigned to \code{Li} as
#' required by the asymmetric projection logic in \code{dist_perpendicular} and \code{dist_parallel}.
#'
#' @inheritParams dist_perpendicular
#' @param w_perp,w_par,w_angle Non-negative weights for the three components. All default to 1.
#' @return Non-negative scalar distance.
#' @export
#' @examples
#' dist_segments(c(0,0), c(10,0), c(2,3), c(8,3))
#' dist_segments(c(0,0), c(10,0), c(0,0), c(10,0))  # identical -> 0
dist_segments <- function(si, ei, sj, ej,
                          w_perp  = 1,
                          w_par   = 1,
                          w_angle = 1) {
  # The paper requires Li to be the longer segment
  len_i <- sum((ei - si)^2)
  len_j <- sum((ej - sj)^2)

  if (len_j > len_i) {
    # Swap so Li is always the longer one
    tmp_s <- si; tmp_e <- ei
    si <- sj;    ei <- ej
    sj <- tmp_s; ej <- tmp_e
  }

  w_perp  * dist_perpendicular(si, ei, sj, ej) +
    w_par   * dist_parallel(si, ei, sj, ej)      +
    w_angle * dist_angle(si, ei, sj, ej)
}
