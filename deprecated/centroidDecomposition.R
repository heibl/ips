## This code is part of the ips package
## © C. Heibl 2016 (2016-11-16)

#' @title Centroid Decomposition of Phylogenetic Trees 
#' @description Splits a phylogenetic tree into a number of subtree
#' with size at most k.
#' @param phy An object of class \code{\link{phylo}}.
#' @param k An integer giving the upper size limit of subtrees.
#' @export

centroidDecomposition <- function(phy, k = 200){
  
  ## initial checks
  ## --------------
  if ( !inherits(phy, "phylo") ) stop("phy is not of class 'phylo'")
  if ( Ntip(phy) <= k ) return(phy) ## nothing to do
  # this is somewhat dirty!
  if ( !is.binary.tree(phy) ) phy <- multi2di(phy, random = FALSE)
  
  n <- Ntip(phy)
  phy <- list(phy)
  
  while ( any(n > k) ) {
    phy <- lapply(phy, splitEqualSubtrees, k = k)
    phy <- unlistFirstLevel(phy)
    phy
    n <- sapply(phy, Ntip)
  }
  phy
}