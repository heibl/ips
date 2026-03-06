## This code is part of the ips package
## © C. Heibl 2014 (last update 2019-11-27)

#' @title Identify/Delete Spurious Rows and Columns from DNA Alignments
#' @description After subsetting (see e.g. \code{\link[ape]{DNAbin}}), DNA sequence
#'   alignments can contain rows and columns that consist entirely of missing
#'   and/or ambiguous character states. \code{identifyEmptyCells} will identify
#'   and \code{deleteEmptyCells} will delete all such rows (taxa) and columns
#'   (characters) from a DNA sequence alignment.
#' @param DNAbin An object of class \code{\link[ape]{DNAbin}}.
#' @param margin A vector giving the subscripts the function will be applied
#'   over: \code{1} indicates rows, \code{2} indicates columns, and \code{c(1,
#'   2)} indicates rows and columns.
#' @param nset A vector of mode character; rows or columns that consist
#'   \strong{only} of the characters given in \code{nset} will be deleted from
#'   the alignment. Allowed are \code{"-"}, \code{"?"},\code{"n"}, \code{"b"},
#'   \code{"d"},\code{"h"}, \code{"v"}, \code{"r"},\code{"y"}, \code{"s"},
#'   \code{"w"},\code{"k"}, and \code{"m"}.
#' @param quiet Logical: if set to \code{TRUE}, screen output will be
#'   suppressed.
#' @details For faster execution, \code{deleteEmptyCells} handles sequences in
#'   \pkg{ape}'s bit-level coding scheme.
#' @return An object of class \code{\link[ape]{DNAbin}}.
#' @references Cornish-Bowden, A. 1984. Nomenclature for incompletely specified
#'   bases in nucleic acid sequences: recommendations 1984. \emph{Nucleic Acids
#'   Res.} \strong{13}: 3021–3030.
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
#'   # identify those rows and colums
#'   (id <- identifyEmptyCells(x))
#'   xx <- x[-id$row, -id$col]
#'   # delete those rows and colums
#'   x <- deleteEmptyCells(x)
#'   image(x)
#'   identical(x, xx)
#' @name EmptyCells


NULL
