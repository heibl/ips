## This code is part of the ips package
## © C. Heibl 2014 (last update 2016-11-16)

#' @title Delete Missing Data from DNA Sequences
#' @description Remove gaps ("-") and/or missing and ambiguous data ("N", "?") 
#' from a sample of DNA sequences.
#' @param x A matrix, a list, or a vector of class \code{\link{DNAbin}}
#' containing the DNA sequences.
#' @return A list or a vector of class \code{DNAbin}.
#' @export

del.miss <- function (x) 
{
  deleteMissing <- function(x) {
    i <- which(x %in% as.raw(c(240, 2, 4)))
    if (length(i)) 
      x[-i]
    else x
  }
  if (!inherits(x, "DNAbin")) 
    x <- as.DNAbin(x)
  if (is.matrix(x)) {
    n <- dim(x)[1]
    y <- vector("list", n)
    for (i in 1:n) y[[i]] <- x[i, ]
    names(y) <- rownames(x)
    x <- y
    rm(y)
  }
  if (!is.list(x)) 
    return(deleteMissing(x))
  x <- lapply(x, deleteMissing)
  class(x) <- "DNAbin"
  x
}