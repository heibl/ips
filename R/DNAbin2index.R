## This code is part of the ips package
## © C. Heibl 2016 (last update 2025-09-14)

#' @title Conversion of DNAbin to Index
#' @description Extract the indices of non-empty positions in a sample of DNA
#'   sequences to zzz.
#' @param x A matrix of class \code{\link[ape]{DNAbin}}.
#' @seealso \code{\link{index2DNAbin}}
#' @export

DNAbin2index <- function(x){
  
  iupac <- c(a = 136, c = 40, g = 72, t = 24, 
             r = 192, y = 48, s = 96, w = 144, k = 80, m = 160, 
             b = 112, d = 208, h = 176, v = 224)
  iupac <- as.raw(iupac)
  
  x <- apply(x, 1, function(msa) which(msa %in% iupac))
  if ( is.matrix(x) ){
    x <- apply(x, 2, as.list)
    x <- lapply(x, unlist)
  }
  x
}