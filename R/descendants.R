## This code is part of the IPS package
## © C. Heibl 2014 (last update 2015-03-26)

descendants <- function(phy, node, type = "t", ignore.tip = TRUE, labels = FALSE){
	
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
  type <- match.arg(type, c("both", "internal", "terminal"))
	tips <- setdiff(edge[, 2], edge[, 1])
  
  # 'node' is a tip 
  # ---------------
  if ( node <= max(tips) ){
    if ( ignore.tip ) x <- node
    else stop("node ", node, " is not an internal node") 
  }
  # normal procedure when 'node' is internal
  # ----------------------------------------
  else {
    x <- edge[edge[,1] == node, 2]
    repeat{
      xx <- x
      x <- sort(unique(c(x, edge[,2][edge[,1] %in% x])))
      if (identical(x, xx)) break
    }
    
    # apply 'type' argument:
    # -----------------------------------------
    if (type == "terminal") {
      x <- intersect(x, tips)
      if (labels) {
        x <- phy$tip.label[x]
      }
    }
    if (type == "internal") x <- setdiff(x, tips)
  }
	return(x)
}