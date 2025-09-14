## PACKAGE: ips
## CALLED BY: mafft
## AUTHOR: Christoph Heibl (at gmx.net)
## LAST UPDATE: 2025-08-30

#' @importFrom ape drop.tip
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