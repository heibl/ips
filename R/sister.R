## This code is part of the megaptera package
## © C. Heibl 2012 (last update 2016-09-21)

sister <- function(phy, node, type = "terminal"){
  
  # checks and definitions
  # ----------------------
  if ( !inherits(phy, "phylo") ) 
    stop ("'phy' is not of class 'phylo'")
  if ( is.character(node) )
    node <- which(phy$tip.label %in% node)
  
  orig.node <- node
  
  if ( length( node) > 1 ){
    if ( !is.monophyletic(phy, node) ){
      stop("elements of 'node' must be monophyletic")
    } else {
      ## Alternative 1: noi (susceptible to duplicate tip labels)
      # node  <- noi(phy, node)
      ## Alternative 2: a small loop ...
      nn <- node
      repeat {
        nn <- sort(unique(phy$edge[phy$edge[, 2] %in% nn, 1]))
        gg <- descendants(phy, min(nn))
        if ( all(node %in% gg) ) break
      } 
      node <- min(nn)
    }
  }
  if ( node == (Ntip(phy) + 1) ){
    stop("node = ", node, " is root node")
    # obj <- descendants(phy, node, type = type)
  } else {
    obj <- phy$edge[, 1][phy$edge[, 2] == node] # getmrca
    obj <- descendants(phy, obj, type = type) # get whole sister clade
    obj <- setdiff(obj, orig.node) # eliminate node 
  }
  obj
}
