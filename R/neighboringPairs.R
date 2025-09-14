## This code is part of the ips package
## Written by C. Heibl 2016 (last update 2025-09-14)

#' @title Neighboring Nodes in a Minimum Spanning Tree
#' @description Finds all pairs of adjacent nodes, i.e. 
#' nodes separated by only one edge, in a minimum spanning 
#' tree
#' @param mst An object of class \code{\link[ape]{mst}}.
#' @export
#' @importFrom ape mst

neighboringPairs <- function(mst){
  
  mst[upper.tri(mst)] <- 0
  id <- which(mst == 1, arr.ind = TRUE)
  data.frame(a = rownames(mst)[id[, 1]],
             b = colnames(mst)[id[, 2]])
}