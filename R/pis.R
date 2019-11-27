## This code is part of the ips package
## © C. Heibl 2014 (last update 2019-11-14)

#' @title Number of Potentially-Informative Sites
#' @description Returns the number or positions of potentially-informative
#'   (parsimony-informative, phylogenetically-informative) sites in DNA sequence
#'   alignment.
#' @param x An object of class \code{\link{DNAbin}}.
#' @param what Either of \code{"absolute"}, \code{"fraction"}, or
#'   \code{"index"}, which will return the absolute number, the relative number
#'   or the indeces of the potentially-informative sites.
#' @param use.ambiguities \emph{Not yet available}.
#' @return Numeric (depending on \code{what}, the number, fraction, or indices
#'   of potentially-informative nucleotide sites).
#' @examples
#' data(ips.16S)	
#' 
#' # number of potentially-informative sites:
#' pis(ips.16S, what = "abs")
#' 
#' # proportion of potentially-informative sites:
#' pis(ips.16S, what = "frac")
#' 
#' # indices of potentially-informative sites:
#' pis(ips.16S, what = "ind")
#' @export

pis  <- function(x, what = "fraction", use.ambiguities = FALSE){

  if (!inherits(x, "DNAbin")) 
    stop("'x' is not of class 'DNAbin'")
  
  what <- match.arg(what, c("absolute", "fraction", "index"))
	
  if (use.ambiguities){
    warning("'use.ambiguities' is currently ignored ", 
            "and IUPAC ambiguity symbols are treated as missing data")
    use.ambiguities <- FALSE
  }
    
  pars.inf <- function(x){
		x <- table(x)
		x <- x[x > 1] # drop apomorphic chars
		n <- c("-", "n", "b", "h", "d", "v", "k", "s", "r", "w", "y")
		if (length(x[!names(x) %in% n]) > 1)				
      x  <-  TRUE									
    else 
      x  <-  FALSE
	}
	x  <-  as.character(x)
	out  <-  apply(x, 2, pars.inf)
  if ( what %in% c("absolute", "fraction") ){
    out <- length(out[out])
    if (what == "fraction"){
      out <- round(out / ncol(x) * 100, digits = 2)
    } 
  } else {
    out <- which(out)
  }
	out
}