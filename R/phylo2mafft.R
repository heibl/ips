## This code is part of the ips package
## © C. Heibl 2014 (last update 2016-11-22)

#' @importFrom utils write.table
#' @export

phylo2mafft <- function(phy){
  
  obj <- matrix(ncol = 2, nrow = 0)
  repeat{
    x <- terminal.clades(phy)
    x <- lapply(x, function(x, phy) phy$tip.label[x], 
                phy = phy)
    x <- lapply(x, sort)
    x <- do.call(rbind, x)
    obj <- rbind(obj, x)
    if ( Ntip(phy) < 3 ) break
    phy <- drop.tip(phy, as.character(x[, 2]))
  }
  obj <- data.frame(obj, 1, 1)
  write.table(obj, "tree.mafft", col.names = FALSE,
              row.names = FALSE, sep = "\t")
  obj
}