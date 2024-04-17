
## This code is part of the ips package
## © C. Heibl 2014 (last update 2016-12-05)

#' @title Identification of Sister Nodes and Clades
#' @description For any given internal node in a phylogeny, this function 
#'   returns the sister clade.
#' @param phy An object of class \code{\link{phylo}}.
#' @param node A vector of mode \code{"numeric"} or \code{"character"} giving 
#'   the number(s) or name(s) of the tiplabel(s); these must be monophyletic.
#' @param type A character string, may be \code{"terminal"}, \code{"internal"}, 
#'   \code{"daughter"}, \code{"all"}, or any unambiguous abbreviation of these; 
#'   \code{"daughter"} will return the MRCA of the sister clade of 
#'   \code{"node"}.
#' @param label Logical, determining if tip number or tip labels will be 
#'   returned.
#' @return A vector of mode \code{"numeric"} or \code{"character"}, containing 
#'   either tip numbers or labels, respectively.
#' @seealso \code{\link{descendants}}, \code{\link{noi}}.
#' @examples 
#' # A phylogeny of bark beetles ...
#' data(ips.tree)
#' tcol <- rep("black", Ntip(ips.tree))
#' tcol[ips.tree$tip.label %in% c("Ips_typographus", "Ips_nitidus")] <- "blue"
#' tcol[ips.tree$tip.label %in% c("Ips_duplicatus")] <- "red"
#' plot(ips.tree, no.margin = TRUE, tip.color = tcol)
#' # What is the sister species of Ips typographus?
#' sister(ips.tree, "Ips_typographus", label = TRUE)
#' # Return the MRCA of the sister clade of Ips duplicatus
#' x <- sister(ips.tree, "Ips_duplicatus", "daughter")
#' nodelabels(node = x, pch = 21, bg = "red")
#' @export

sister <- function(phy, node, type = "terminal", label = FALSE){
  
  # checks and definitions
  # ----------------------
  if ( !inherits(phy, "phylo") ) 
    stop ("'phy' is not of class 'phylo'")

  if ( is.character(node) | is.factor(node) )
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
  if ( label ) obj <- phy$tip.label[obj]
  obj
}
