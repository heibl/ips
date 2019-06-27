## This code is part of the ips package
## © C. Heibl 2014 (last update 2019-06-20)

#' @title Equal Tip Heights
#' @description Modify terminal edge lengths to create "exactly" (see Details)
#'   equal tip heights (sum of edge lengths from root to tip)
#' @param phy An object of class \code{\link[ape]{phylo}}.
#' @param baseline A character string giving a function to calculate the
#'   baseline tip height, e.g. \code{"min"}, \code{"max"} or \code{"mean"}.
#' @details What is "exactly" equal depends on the precision of the system
#'   (\code{\link{.Machine}}); in any case the resulting phylogeny will pass
#'   \code{\link{is.ultrametric}} with default arguments.
#' @return An object of class \code{\link[ape]{phylo}} with changed terminal
#'   edge lengths.
#' @note \code{forceEqualTipHeights} is only intended to correct small rounding
#'   errors in edge lengths, not to make an additive phylogeny ultrametric. For
#'   the latter, see e.g. \code{\link{chronos}}.
#' @seealso \code{\link{tipHeights}}
#' @export

forceEqualTipHeights <- function(phy, baseline = "mean"){
  
  tip_heights <- tipHeights(phy)
  message("Range of tip heights is ", 
      format((max(tip_heights) - min(tip_heights)) / max(tip_heights) * 100, scientific = FALSE),
      "% of maximum tip height")
  baseline <- eval(parse(text = paste0(baseline, "(tip_heights)")))
  message("\nBaseline tip height: ", baseline)
  diffs <- tip_heights - baseline
  diffs <- data.frame(diff = diffs[diffs != 0])
  if (!nrow(diffs)) stop("tip heights are already equal")
  diffs$tip_edge <- which(phy$edge[, 2] %in% which(phy$tip.label %in% rownames(diffs)))
  for (i in 1:nrow(diffs)){
    phy$edge.length[diffs$tip_edge[i]] <- phy$edge.length[diffs$tip_edge[i]] - diffs$diff[i]
  }
  phy
}


