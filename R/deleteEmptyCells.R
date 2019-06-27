## This code is part of the ips package
## © C. Heibl 2014 (last update 2019-06-20)

#' @title Delete Spurious Rows and Columns from DNA Alignments
#' @description After subsetting (see e.g. \code{\link{DNAbin}}), DNA sequence
#'   alignments can contain rows and columns that consist entirely of missing
#'   and/or ambiguous character states. \code{deleteEmptyCells} will delete all
#'   such rows (taxa) and columns (characters) from a DNA sequence alignment.
#' @param DNAbin An object of class \code{\link{DNAbin}}.
#' @param margin A vector giving the subscripts the function will be applied
#'   over: \code{1} indicates rows, \code{2} indicates columns, and
#'   \code{c(1, 2)} indicates rows and columns.
#' @param nset A vector of mode character; rows or columns that consist 
#'   \strong{only} of the characters given in \code{nset} will be deleted from 
#'   the alignment. Allowed are \code{"-"}, \code{"?"},\code{"n"}, \code{"b"}, 
#'   \code{"d"},\code{"h"}, \code{"v"}, \code{"r"},\code{"y"}, \code{"s"}, 
#'   \code{"w"},\code{"k"}, and \code{"m"}.
#' @param quiet Logical: if set to \code{TRUE}, screen output will be suppressed.
#' @details For faster execution, \code{deleteEmptyCells} handles sequences in 
#'   \pkg{ape}'s bit-level coding scheme. 
#' @return An object of class \code{\link{DNAbin}}.
#' @references Cornish-Bowden, A. 1984. Nomenclature for incompletely 
#' specified bases in nucleic acid sequences: recommendations 1984. 
#' \emph{Nucleic Acids Res.} \strong{13}: 3021–3030.
#' @seealso \code{\link{trimEnds}}, \code{\link{deleteGaps}}
#' @examples 
#'   # COX1 sequences of bark beetles
#'   data(ips.cox1)
#'   # introduce completely ambiguous rows and colums
#'   x <- as.character(ips.cox1[1:6, 1:60])
#'   x[3, ] <- rep("n", 60)
#'   x[, 20:24] <- rep("-", 6)
#'   x <- as.DNAbin(x)
#'   image(x)
#'   # delete those rows and colums
#'   x <- deleteEmptyCells(x)
#'   image(x)
#' @export

deleteEmptyCells <- function(DNAbin, margin = c(1, 2),
                             nset = c("-", "n", "?"),
                             quiet = FALSE){
  
  if ( !inherits(DNAbin, "DNAbin") ) 
    stop("'DNAbin' is not of class 'DNAbin'")
  
  ## convert character to raw
  
  ## IUPAC ambiguity code
  ## --------------------
  iupac <- c(n = 240, "?" = 2, "-" = 4,
    # a = 136, c = 40, g = 72, t = 24, 
    r = 192, y = 48, s = 96, w = 144, k = 80, m = 160, 
    b = 112, d = 208, h = 176, v = 224)
  nset <- iupac[nset]
  nset <- as.raw(nset)
  
  ## function that detects non-empty strings
  isNotEmpty <- function(x, nset){
    ifelse(all(unique(x) %in% nset), FALSE, TRUE)
  }
  
  size <- dim(DNAbin)
  
  ## rows  (margin == 1)
  if (1 %in% margin){
    rowind <- which(apply(DNAbin, 1, isNotEmpty, nset = nset))
    DNAbin <- DNAbin[rowind, ]
  }
  
  ## columns (margin == 2)
  if (2 %in% margin){
    colind <- which(apply(DNAbin, 2, isNotEmpty, nset = nset))
    DNAbin <- DNAbin[, colind]
  }
  
  ## screen output (if desired)
  if (!quiet) {
    size <- size - dim(DNAbin)
    rows <- ifelse(size[1] == 1, " row ", " rows ")
    cols <- ifelse(size[2] == 1, " column ", " columns ")
    message(size[1], rows, "deleted from alignment\n",
            size[2], cols, "deleted from alignment")
  }  
  DNAbin
}

