## This code is part of the ips package
## © C. Heibl 2016 (last update 2016-11-23)

#' @title Numbers of Tips of (Sub)trees
#' @description Counts the number of tips of a given clade of a 
#' phylogenetic tree.
#' @param phy An object of class \code{\link{phylo}}.
#' @param node An integer given the number of an internal node.
#' @return An integer giving the number of tips.
#' @examples 
#' set.seed(1234)
#' tr <- rtree(12)
#' plot(tr); nodelabels()
#' ntip(tr, 16)
#' @export

ntip <- function(phy, node){
  
  nmax <- Ntip(phy)
  tips <- vector()
  while ( length(node) > 0 ){
    node <- phy$edge[phy$edge[, 1] %in% node, 2]
    tips <- c(tips, node[node <= nmax])
    node <- node[node > nmax]
  }
  length(tips)
}