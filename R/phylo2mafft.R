## This code is part of the ips package
## © C. Heibl 2014 (last update 2019-07-04)

#' @title Convert Trees for MAFFT
#' @description Converts a phylogenetic tree of class \code{"phylo"} to a format
#'   usable as a guide tree by MAFFT. This function is called internally by 
#'   \code{\link{mafft}}.
#' @param phy A phylogenetic tree of class \code{\link{phylo}}.
#' @param file A character string giving a filename. May be missing, in which
#'   case the results are only printed on the screen.
#' @return A matrix coding the MAFFT-formatted tree, as a side effect the same 
#'   matrix is written to \code{file}.
#' @references The MAFFT website:
#'   \url{https://mafft.cbrc.jp/alignment/software/index.html}
#' @seealso \code{\link{mafft}} for an interface to MAFFT.
#' @importFrom ape is.binary multi2di
#' @importFrom utils write.table
#' @export

phylo2mafft <- function(phy, file){
  
  if (!is.binary(phy)) phy <- multi2di(phy)

  obj <- matrix(ncol = 2, nrow = 0)
  repeat {
    x <- terminal.clades(phy)
    x <- lapply(x, function(x, phy) phy$tip.label[x],
                phy = phy)
    x <- lapply(x, sort)
    x <- do.call(rbind, x)
    obj <- rbind(obj, x)
    if (Ntip(phy) < 3) break
    phy <- drop.tip(phy, as.character(x[, 2]))
  }
  obj <- data.frame(obj, 1, 1)
  if (!missing(file)){
    write.table(obj, file, col.names = FALSE,
                row.names = FALSE, sep = "\t")
  }
  obj
}
