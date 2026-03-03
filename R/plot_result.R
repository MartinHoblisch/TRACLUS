#' @importFrom graphics segments legend lines points
NULL

#' Plot TRACLUS clustering results
#'
#' Visualises the output of \code{cluster_segments} or \code{traclus}: noise segments in grey,
#' clustered segments and their representative trajectories in colour. Two display modes are
#' available (see \code{mode}).
#'
#' @param clustered_segs Data frame with columns \code{sx}, \code{sy}, \code{ex}, \code{ey},
#'   and \code{cluster_id} (as returned by \code{cluster_segments} or \code{traclus()$segments}).
#' @param representatives Named list of representative trajectory matrices
#'   (as returned by \code{compute_all_representatives} or \code{traclus()$representatives}).
#'   If \code{NULL}, only segments are plotted.
#' @param mode Character string controlling how clusters are displayed:
#'   \describe{
#'     \item{\code{"clusters"}}{Each cluster's segments and representative trajectory share the
#'       same colour. The representative is drawn with a black outline to stand out from
#'       the underlying segments. This is the default.}
#'     \item{\code{"representatives"}}{Clustered segments are drawn in a uniform grey (solid,
#'       not dashed). Only the representative trajectories are coloured, making them the
#'       visual focus.}
#'   }
#' @param title Plot title. Default \code{"TRACLUS clustering result"}.
#' @param noise_col Colour for noise segments. Default \code{"grey80"}.
#' @param cluster_col Colour for clustered segments in \code{"representatives"} mode.
#'   Default \code{"grey60"}.
#' @param noise_lwd Line width for noise segments. Default 0.5.
#' @param seg_lwd Line width for clustered segments. Default 1.5.
#' @param rep_lwd Line width for representative trajectories. Default 4.
#' @param pal Character vector of colours used for clusters / representatives. Cycles if there
#'   are more clusters than colours. Defaults to the Wong (2011) colorblind-safe palette.
#' @param show_legend Logical. If \code{TRUE} (default), a legend is drawn.
#' @return Invisibly \code{NULL}, called for its side effect.
#' @export
#' @examples
#' result <- traclus(make_toy_trajectories(), eps = 10, min_lns = 3)
#' plot_traclus_result(result$segments, result$representatives)
#' plot_traclus_result(result$segments, result$representatives, mode = "representatives")
plot_traclus_result <- function(clustered_segs,
                                representatives = NULL,
                                mode        = c("clusters", "representatives"),
                                title       = "TRACLUS clustering result",
                                noise_col   = "grey80",
                                cluster_col = "grey60",
                                noise_lwd   = 0.5,
                                seg_lwd     = 1.5,
                                rep_lwd     = 4,
                                pal         = NULL,
                                show_legend = TRUE) {
  mode <- match.arg(mode)

  stopifnot(
    is.data.frame(clustered_segs),
    all(c("sx", "sy", "ex", "ey", "cluster_id") %in% names(clustered_segs))
  )

  if (is.null(pal)) {
    pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
             "#0072B2", "#D55E00", "#CC79A7", "#000000")
  }

  # Coordinate range across all segments
  all_x <- c(clustered_segs$sx, clustered_segs$ex)
  all_y <- c(clustered_segs$sy, clustered_segs$ey)

  plot(NULL,
       xlim = range(all_x), ylim = range(all_y),
       xlab = "x (km)", ylab = "y (km)",
       main = title, asp = 1)

  is_noise    <- is.na(clustered_segs$cluster_id)
  cluster_ids <- sort(unique(clustered_segs$cluster_id[!is_noise]))

  # --- Noise segments (always dashed grey) ---
  if (any(is_noise)) {
    ns <- clustered_segs[is_noise, ]
    segments(ns$sx, ns$sy, ns$ex, ns$ey,
             col = noise_col, lwd = noise_lwd, lty = 2)
  }

  # --- Clustered segments ---
  if (mode == "clusters") {
    # Each cluster in its own colour
    for (k in seq_along(cluster_ids)) {
      cid  <- cluster_ids[k]
      csub <- clustered_segs[!is_noise & clustered_segs$cluster_id == cid, ]
      col_k <- pal[(k - 1L) %% length(pal) + 1L]
      segments(csub$sx, csub$sy, csub$ex, csub$ey,
               col = col_k, lwd = seg_lwd)
    }
  } else {
    # All clustered segments in uniform grey (solid, not dashed)
    cs <- clustered_segs[!is_noise, ]
    if (nrow(cs) > 0L) {
      segments(cs$sx, cs$sy, cs$ex, cs$ey,
               col = cluster_col, lwd = seg_lwd)
    }
  }

  # --- Representative trajectories ---
  if (!is.null(representatives) && length(representatives) > 0L) {
    for (k in seq_along(representatives)) {
      rt    <- representatives[[k]]
      col_k <- pal[(k - 1L) %% length(pal) + 1L]

      if (mode == "clusters") {
        # Black outline: draw a wider black line first, then the coloured line on top
        lines(rt[, "x"], rt[, "y"], col = "black", lwd = rep_lwd + 2)
        lines(rt[, "x"], rt[, "y"], col = col_k,   lwd = rep_lwd)
        points(rt[, "x"], rt[, "y"], col = col_k, pch = 19, cex = 0.6)
      } else {
        lines(rt[, "x"], rt[, "y"], col = "black", lwd = rep_lwd + 2)
        lines(rt[, "x"], rt[, "y"], col = col_k,   lwd = rep_lwd)
        points(rt[, "x"], rt[, "y"], col = col_k, pch = 19, cex = 0.6)
      }
    }
  }

  # --- Legend ---
  if (show_legend && length(cluster_ids) > 0L) {
    legend_labels <- paste("Cluster", cluster_ids)
    legend_cols   <- pal[(seq_along(cluster_ids) - 1L) %% length(pal) + 1L]
    legend("bottomright",
           legend = legend_labels,
           col    = legend_cols,
           lwd    = 3, cex = 0.8, bg = "white")
  }

  invisible(NULL)
}
