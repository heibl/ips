## This code is part of the ips package
## © C. Heibl 2014 (last update 2017-11-28)

#' @export

tipHeights <- function(phy){
  
  if (!inherits(phy, "phylo")) stop("'phy' is not of class 'phylo'")
  if (is.null(phy$edge.length)) stop("'phy' has no branch lengths")
  
  ## sum up heights from tip to root
  ## -------------------------------
  nodes <- 1:Ntip(phy)
  heights <- rep(0, Ntip(phy))
  while ( !all(is.na(nodes)) ){
    id <- match(nodes, phy$edge[, 2])
    h <- phy$edge.length[id]
    h[is.na(h)] <- 0
    heights <- heights + h
    nodes <- phy$edge[id, 1]
  }
  
  ## name vector elements (works also with unorthodox tip numbering)
  tips <- phy$edge[phy$edge[, 2] %in%  1:Ntip(phy), 2]
  names(heights) <- phy$tip.label[tips]
  
  heights
}