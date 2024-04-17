
## This code is part of the ips package
## © C. Heibl 2014 (last update 2016-11-07)

#' @export

descendants <- function(phy, node, type = "t", ignore.tip = TRUE, 
                        labels = FALSE){
	
  # checks and definitions
  # ----------------------
  if ( inherits(phy, "phylo") ){
    edge <- phy$edge
  } else {
    if ( !is.matrix(phy) ) {
      stop("'phy' must be of classes 'phylo' or 'matrix'")
    } else {
      edge <- phy
      labels <- FALSE
    }
  }
  if ( length(node) > 1) stop("'node' must be vector of length 1")
  type <- match.arg(type, c("all", "daughter", "internal", "terminal"))
	tips <- setdiff(edge[, 2], edge[, 1])
  
  # 'node' is a tip 
  # ---------------
  if ( node <= max(tips) ){
    if ( ignore.tip ){
      x <- node
    } else {
      stop("node ", node, " is not an internal node") 
    }
  } else {
    
    # normal procedure when 'node' is internal
    # ----------------------------------------
    x <- edge[edge[,1] == node, 2] # immediate daughter nodes
    if ( type %in% c("internal", "terminal", "all") ){
      repeat{
        xx <- x
        x <- sort(unique(c(x, edge[,2][edge[,1] %in% x])))
        if (identical(x, xx)) break
      }
      if ( type == "internal" ) x <- setdiff(x, tips)
    }
  }
	## apply 'type' argument:
	## ----------------------
	if ( type == "terminal" ) {
	  x <- intersect(x, tips)
	  if (labels) {
	    x <- phy$tip.label[x]
	  }
	}
	x
}