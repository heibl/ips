## This code is part of the ips package
## © C. Heibl 2014 (last update 2016-11-29)

#' @title Find Monophyletic Subsets in Species Lists
#' @description Takes a phylogeny and a subset of its tiplabels
#' and splits the list of tiplabels into monophyletic groups 
#' (clades).
#' @param phy An object of class \code{\link{phylo}}.
#' @param tips A vector of mode \code{"character"} containing any 
#' subset of the tiplabels in \code{phy}.
#' @return A list.
#' @export

splitIntoClades <- function(phy, tips){
  
  obj <- vector()
  repeat {
    ss <- lapply(tips, sister, phy = phy, label = TRUE)
    id <- sapply(ss, function(z, s) all(z %in% s), s = tips)
    
    obj <- c(obj, tips[!id])
    if ( !any(id) ) break
    ss <- ss[id]
    tips <- tips[id]
    
    for ( i in seq_along(tips) ) ss[[i]] <- sort(c(ss[[i]], tips[i]))
    tips <- unique(ss)
    
    ## check if any one element is a subset of another
    id <- vector(length = length(tips))
    for ( i in seq_along(tips) ){
      z <- lapply(tips[-i], setdiff, x = tips[[i]])
      id[i] <- !any(sapply(z, length) == 0)
    }
    tips <- tips[id]
  }
  obj
}