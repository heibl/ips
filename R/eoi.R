## This code is part of the ips package
## © C. Heibl 2014 (last update 2017-12-20)

#' @rdname oi
#' @export

eoi <- function(phy, node, group, regex = FALSE, stem = FALSE, 
                monophyletic = FALSE){
  
  ## check phylogeny
  if ( !inherits(phy, "phylo") ) 
    stop("'phy' is not of class 'phylo'")  
  
  ## if group argument is specified, get node numbers first
  ## ------------------------------------------------------
  if (missing(node)){
    if (missing(group)) stop("either 'node' or 'group' must be specified")
    node <- noi(phy = phy, group = group, regex = regex, stem = stem, 
                monophyletic = monophyletic)
  }
  
  ## get edge numbers from node numbers
  ## ----------------------------------
  sapply(node, function(phy, z) which(phy$edge[, 2] == z), phy = phy)
}