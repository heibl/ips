## This code is part of the ips package
## © C. Heibl 2014 (last update 2019-06-27)

#' @title Masking of Sequence Alignments with GBLOCKS
#' @description Provides a wrapper to Gblocks, a computer program written in
#'   ANSI C language that eliminates poorly aligned positions and divergent
#'   regions of an alignment of DNA or protein sequences. Gblocks selects
#'   conserved blocks from a multiple alignment according to a set of features
#'   of the alignment positions.
#' @param x A matrix of DNA sequences of classes \code{\link[ape]{DNAbin}}.
#' @param b1 A real number, the \bold{minimum number of sequences for a
#'   conserved position} given as a fraction. Values between 0.5 and 1.0 are
#'   allowed. \emph{Larger} values will \emph{decrease} the number of selected
#'   positions, i.e. are more \emph{conservative}. Defaults to 0.5
#' @param b2 A real number, the \bold{minimum number of sequences for a flank
#'   position} given as a fraction. Values must be equal or larger than
#'   \code{b1}. \emph{Larger} values will \emph{decrease} the number of selected
#'   positions, i.e. are \emph{more conservative}. Defaults to 0.5
#' @param b3 An integer, the \bold{maximum number of contiguous nonconserved
#'   positions}; any integer is allowed. \emph{Larger} values will
#'   \emph{increase} the number of selected position, i.e. are \emph{less
#'   conservative}. Defaults to the number of positions in the alignment.
#' @param b4 An integer, the \bold{minimum length of a block}, any integer equal
#'   to or bigger than 2 is allowed. \emph{Larger} values will \emph{decrease}
#'   the number of selected positions, i.e. are \emph{more conservative}. Defaults to
#'   2.
#' @param b5 A character string indicating the \bold{treatment of gap
#'   positions}. Three choices are possible. 1. \code{"n"}: \emph{No} gap
#'   positions are allowed in the final alignment. All positions with a single
#'   gap or more are treated as a gap position for the block selection
#'   procedure, and they and the adjacent nonconserved positions are eliminated.
#'   2. \code{"h"}: Only positions where \emph{50\% or more} of the sequences
#'   have a gap are treated as a gap position. Thus, positions with a gap in
#'   less than 50\% of the sequences can be selected in the final alignment if
#'   they are within an appropriate block. 3. \code{"a"}: \emph{All} gap
#'   positions can be selected. Positions with gaps are not treated differently
#'   from other positions (default).
#' @param target A vector of mode \code{"character"} giving the output format:
#'   \code{"alignment"} will return the alignment with only the selected
#'   positions, \code{"index"} will return the indices of the selected position,
#'   and \code{"score"} will provide a score for every position in the original
#'   alignment (0 for excluded, 1 for included).
#' @param exec A character string indicating the path to the GBLOCKS executable.
#' @details Explanation of the routine taken from the Online Documentation:
#'   First, the degree of conservation of every positions of the multiple
#'   alignment is evaluated and classified as \emph{nonconserved},
#'   \emph{conserved}, or \emph{highly conserved}. All stretches of contiguous
#'   nonconserved positions bigger than a certain value (\bold{b3}) are
#'   rejected. In such stretches, alignments are normally ambiguous and, even
#'   when in some cases a unique alignment could be given, multiple hidden
#'   substitutions make them inadequate for phylogenetic analysis. In the
#'   remaining blocks, flanks are examined and positions are removed until
#'   blocks are surrounded by highly conserved positions at both flanks. This
#'   way, selected blocks are anchored by positions that can be aligned with
#'   high confidence. Then, all gap positions -that can be defined in three
#'   different ways (\bold{b5})- are removed. Furthermore, nonconserved
#'   positions adjacent to a gap position are also eliminated until a conserved
#'   position is reached, because regions adjacent to a gap are the most
#'   difficult to align. Finally, small blocks (falling below the limit of
#'   \bold{b4}) remaining after gap cleaning are also removed.
#' @return A \code{matrix} of class \code{"DNAbin"}
#' @note \code{gblocks} was last updated and tested to work with Gblocks 0.91b.
#'   If you have problems getting the function to work with a newer version of
#'   Gblocks, please contact the package maintainer.
#' @references Castresana, J. 2000. Selection of conserved blocks from multiple
#'   alignments for their use in phylogenetic analysis. \emph{Molecular Biology
#'   and Evolution} \bold{17}, 540-552.
#' @references Talavera, G., and J. Castresana. 2007. Improvement of phylogenies
#'   after removing divergent and ambiguously aligned blocks from protein
#'   sequence alignments. \emph{Systematic Biology} \bold{56}, 564-577.
#' @references  \bold{Gblocks website}: 
#'   \url{https://www.biologiaevolutiva.org/jcastresana/Gblocks.html}
#' @seealso \code{\link{mafft}} and \code{\link{prank}} for multiple sequence
#'   alignment; \code{\link{aliscore}} for another alignment masking algorithm.
#' @examples
#' data(ips.28S)
#' \dontrun{gblocks(ips.28S)}  
#' @export

gblocks <- function(x, b1 = .5, b2 = b1, b3 = ncol(x), 
                    b4 = 2, b5 = "a", target = "alignment", exec){
  
  ## Check path to executable
  ## ------------------------
  if (missing(exec)){
    ## Try to guess location of executable
    exec <- list.files(path = "/Applications", pattern = "Gblocks")
    exec <- list.files(path = file.path("/Applications", exec), 
                       pattern = "Gblocks", full.names = TRUE)
    if (!length(exec)) stop("path to executable not given")
    message("Using '", exec, "'", appendLF = TRUE)
  } else {
    if (!file.exists(exec)) 
      stop("executable '", exec, "' does not exist", sep = "")
  }
  
  if (inherits(x, "alignment")) stop("cannot handle class 'alignment'")
  if (inherits(x, "list")) stop("cannot handle unaligned sequences")
  
  target <- match.arg(target, c("alignment", "index","score"))
  
  ## check parameters:
  ## -----------------
  if (b1 < .5 | b1 > 1) stop ("b1 not in [0.5, 1]")
  if (b2 < b1 | b2 > 1) stop ("b2 not in [b1, 1]")
  if (b3 < 0 | b4 > ncol(x)) stop ("b3 not in [0, ", ncol(x), "]")
  if (b4 < 2 | b4 > ncol(x)) stop ("b4 not in [2, ", ncol(x), "]")
  b5 <- match.arg(b5, c("a", "h", "n"))
  
  ntax <- nrow(x)
  b1 <- floor(ntax * b1) + 1
  b2 <- floor(ntax * b2) + 1
  
  message("Gblocks parameters:",
          "\n- minimum number of sequences for a conserved position : ", b1,
          "\n- minimum number of sequences for a flank position     : ", b2,
          "\n- maximum number of contiguous nonconserved positions  : ", b3,
          "\n- minimum length of a block                            : ", b4,
          "\n- allowed gap positions                                : ", b5)
  
  write.fas(x, "R2GBLOCK.fas")
  system(paste0(exec, " R2GBLOCK.fas -t=d", 
               " -b1=", b1,
               " -b2=", b2,
               " -b3=", b3,
               " -b4=", b4,
               " -b5=", b5))
  
  if (target == "alignment"){
    out <- read.FASTA("R2GBLOCK.fas-gb")
  } else {
    out <- scan("R2GBLOCK.fas-gb.htm", what = "c", sep = "\n", quiet = TRUE)
    out <- out[grep("^Flanks: ", out)]
    out <- gsub("^Flanks: [[]|[]]  $", "", out)
    out <- gsub("]  [", ",", out, fixed = TRUE)
    out <- gsub("  ", ":", out, fixed = TRUE)
    out <- eval(parse(text = paste0("c(", out, ")")))
    if (target == "score"){
      out <- as.numeric(1:ncol(x) %in% out)
    }
  }
  unlink(list.files(pattern = "R2GBLOCK"))
  out
}