## This code is part of the ips package
## © C. Heibl 2014 (last update 2025-09-14)

#' @title Collapse Edges in a Phylogeny
#' @description Given a set of node support values (e.g., bootstrap proportions,
#'   posterior probabilities) and a certain threshold, all edges receiving less
#'   support than the threshold will be collapsed.
#' @param phy An object of class \code{\link[ape]{phylo}}.
#' @param value A character string giving the name of the list element of the
#'   \code{\link[ape]{phylo}} object that contains the support values; default is
#'   \code{"node.label"}.
#' @param cutoff A numeric value giving the threshold below which edges will be
#'   collapsed.
#' @details Collapsing weekly supported edges can help draw attention to those
#' parts of a phylogeny that are really supported by the underlying data.
#' @return An object of class \code{\link[ape]{phylo}}.
#' @examples
#' ## phylogeny of bark beetles
#' data(ips.tree)
#'
#' ## non-parametric bootstrap proportions (BP)
#' ips.tree$node.label
#'
#' ## collapse clades with < 70 BP
#' tr <- collapseUnsupportedEdges(ips.tree, "node.label", 70)
#'
#' ## show new topology
#' plot(tr, no.margin = TRUE)
#'
#' @importFrom ape Ntip
#' @export



collapseUnsupportedEdges <- function(phy, value, cutoff){
  
  if ( !inherits(phy, "phylo") ) 
    stop ("'phy' is not of class 'phylo'")
  
  if ( missing(value) ) value <- "node.label"
  
  stat <- as.numeric(phy[[value]])
  nt <- Ntip(phy)
  root.node <- nt + 1
  collapse <- which(stat < cutoff) + nt
  ## the root node cannot be collapsed:
  collapse <- setdiff(collapse, root.node)
  
  ## collapse nodes in post-order traversal!!
  ## ----------------------------------------
  for ( i in rev(collapse) ){
    
    #i <- rev(collapse)[1] # FOR DEBUGGING
    
    ## identify edges
    id <- phy$edge[, 2] == i
    id2 <- phy$edge[, 1] == i
    
    ## modify: node values
    #phy$node.label <- phy$node.label[!id]
    
    ## modify: edges
    phy$edge[id2, 1] <- phy$edge[id, 1]
    phy$edge <- phy$edge[!id, ]
    phy$edge[phy$edge > i] <- phy$edge[phy$edge > i] - 1
    
    ## modify: edge lengths
    phy$edge.length[id2] <- phy$edge.length[id2] + phy$edge.length[id]
    phy$edge.length <- phy$edge.length[!id]
    
    ## modify: node labels
    phy[["node.label"]] <- phy[["node.label"]][-(i - nt)]
    
    ## modify: number of internal nodes
    phy$Nnode <- phy$Nnode - 1
    
  }
  phy
}