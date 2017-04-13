## This code is part of the ips package
## © C. Heibl 2014 (last update 2017-03-22)

#' @title Find Pairs of Sister Species
#' @description Finds pairs of sister species in
#' a phylogenetic tree.
#' @param phy An object of class \code{\link{phylo}}.
#' @param labels Logical, indicating whether to return tip labels or tip numbers.
#' @return A list of which each element contains the tip labels
#' of a sister species pair.
#' @examples
#' set.seed(1234)
#' tr <- rtree(12)
#' plot(tr)
#' terminalSisters(tr)
#' @export
#' @import ape

terminalSisters <- function(phy, labels = TRUE){

  obj <- lapply(1:Ntip(phy), sister, phy = phy)
  for ( i in seq_along(obj))
    obj[[i]] <- sort(c(obj[[i]], i))
  obj <- unique(obj)
  is.nested <- function(x, y){
    identical(sort(union(x, y)), x)
  }
  id <- vector(length = length(obj))
  for ( i in seq_along(obj) ){
    id[i] <- !any(sapply(obj[-i], is.nested, x = obj[[i]]))
  }
  obj <- obj[id]
  if (labels) obj <- lapply(obj, function(phy, x) phy$tip.label[x], phy = phy)
  obj
}
