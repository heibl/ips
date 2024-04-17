## This code is part of the ips package
## © C. Heibl 2016 (last update 2019-03-07)

## WORK UNDER PROGRESS
## - cannot handle polytomies
## - does not incoporate branch length information

#' @title Conversion from PHYLO to MST Object
#' @description Converts a phylogenetic tree (class \code{phylo}) into a
#'   minimum spanning tree (class \code{mst}).
#' @details The current version of \code{phylo2mst} does not handle polytomies
#'   and does not incorporate branch length information. Note that topological 
#'   information is lost during the conversion.
#' @param phy An object of class \code{\link{phylo}}.
#' @examples
#' phy <- rtree(12)
#' plot(phy)
#' mst <- phylo2mst(phy)
#' plot(mst)
#' @importFrom ape drop.tip is.binary Ntip
#' @export 

phylo2mst <- function(phy){
  
  ## do some checks:
  if (!is.binary(phy)) stop("cannot handle polytomies")
  
  ## create empty object of class mst:
  mst <- matrix(data = 0, nrow = Ntip(phy), ncol = Ntip(phy),
                dimnames = list(phy$tip.label, phy$tip.label))
  class(mst) <- "mst"
  
  ## doe the conversion:
  breakafter <- FALSE
  repeat {
    tp <- terminalSisters(phy)
    for (j in seq_along(tp)){
      
      if (Ntip(phy) == 2) breakafter <- TRUE
      tpp <- tp[[j]]
      mst[tpp[1], tpp[2]] <- 1
      mst[tpp[2], tpp[1]] <- 1
      if (Ntip(phy) > 2) phy <- drop.tip(phy, tpp[2])
    }
    if (breakafter) break
  }
  
  mst
}