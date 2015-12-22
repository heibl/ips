## This code is part of the ips package
## © C. Heibl 2014 (last update 2015-03-26)

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