## This code is part of the ips package
## © C. Heibl 2014 (last update 2017-12-19)

#' @export

eoi <- function(noi){
  
  ## check phylogeny
  if ( !inherits(phy, "phylo") ) 
    stop("'phy' is not of class 'phylo'")  
  
  sapply(noi, function(phy, z) which(phy$edge[, 2] == z), phy = phy)
}