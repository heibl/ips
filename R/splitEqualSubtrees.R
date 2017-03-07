## This code is part of the ips package
## © C. Heibl 2016 (last update 2016-11-23)

#' @title Split Phylogenetic Trees
#' @description Splits a phylogenetic tree into roughly equally sized
#' subtrees
#' @param phy An object of class \code{\link{phylo}}.
#' @param k An integer, only trees with > k tips will be split.
#' @return A list of two phylogenetic trees.
#' @seealso \code{\link{extract.clade}}.
#' @examples 
#' set.seed(1234)
#' tr <- rtree(12000)
#' subtr <- splitEqualSubtrees(tr)
#' subtr
#' @export

splitEqualSubtrees <- function(phy, k){
  
  ## only split if size is bigger than k
  if ( !missing(k) ){
    if ( Ntip(phy) <= k ) return(phy)
  }
  
  
  repeat {
    
    ## calculate number of tips of root daughters
    ## and their difference
    d <- phy$edge[phy$edge[, 1] == Ntip(phy) + 1, 2]
    n <- sapply(d, ntip, phy = phy)
    if ( n[1] == n[2] ) break
    df <- diff(sort(n))
    
    ## calulate number of tips of descendants of
    ## root descendant with more tips
    dd <- phy$edge[phy$edge[, 1] == d[which.max(n)], 2]
    nn <- sapply(dd, ntip, phy = phy)
    ## root on the granddaughter that has more tips
    phy.new <- root(phy, node = dd[which.max(nn)], resolve.root = TRUE)
    
    ## calculate difference in tips
    d.new <- phy.new$edge[phy.new$edge[, 1] == Ntip(phy) + 1, 2]
    n.new <- sapply(d.new, ntip, phy = phy.new)
    df.new <- diff(sort(n.new))
    
    if ( df.new < df ){
      
      ## accept changes and make new suggestion
      phy <- phy.new
      
    } else {
      
      ## refuse proposal and take previous root position
      break
    }
  }
  
  ## return list with equally sized subtrees
  list(extract.clade(phy, d[1]),
       extract.clade(phy, d[2]))
}
