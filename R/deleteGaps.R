## This code is part of the ips package
## © C. Heibl 2014 (last update 2016-11-16)

#' @title Remove Gap Positions From DNA Sequences
#' @description Remove indel positions (or gaps) from a DNA sequence 
#' alignment. For faster execution, \code{deleteGaps} handles sequences 
#' in \pkg{ape}'s bit-level coding scheme.
#' @param x An object of class \code{\link{DNAbin}}.
#' @param gap.max An integer, which gives the maximum number of
#' gap characters ("-") that will be tolerated at any given alignment 
#' position (column). Only values between \code{0} and \code{nrow(x) - 4}
#' make sense phylogenetically.
#' @details The default, \code{nmax = nrow(x) - 4}, removes all those 
#' positions from the alignment, which contain at least four non-gap 
#' characters, which is the minimum number of sequences needed to 
#' produce a non-trivial unrooted topology. All gaps will be excluded 
#' by selecting \code{nmax = 0} and half of all gaps with \code{nmax = 
#' nrow(x) / 2}.  
#'
#' In contrast, \code{\link[ape]{del.gaps}} removes all gap characters 
#' from the alignment, so most probably the result will not be a set of 
#' sequences of equal length and the matrix will be coerced to a list.
#' @return An object of class \code{\link{DNAbin}.}
#' @seealso  \code{\link{code.simple.gaps}} for coding of simple gaps, 
#' \code{\link[ape]{del.gaps}} for removal of all gap symbols from an 
#' alignment, \code{\link{gblocks}} and \code{\link{aliscore}} for more 
#' sophisticated  methods of cleaning/masking alignments.
#' @export

deleteGaps <- function(x, gap.max = nrow(x) - 4){
  
  if (!inherits(x, "DNAbin")) stop("'x' is not of class 'DNAbin'")
  
  id <- apply(x, 2, function(z) sum(z == 4))
  x[, which(id < gap.max)]
}