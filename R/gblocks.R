## This code is part of the ips package
## © C. Heibl 2014 (last update 2014-11-27)

#' Masking of Sequence Alignments with GBLOCKS
#' 
#' This function provides a wrapper to Gblocks, a computer program written in
#' ANSI C language that eliminates poorly aligned positions and divergent
#' regions of an alignment of DNA or protein sequences. Gblocks selects
#' conserved blocks from a multiple alignment according to a set of features of
#' the alignment positions.
#' 
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
#' @param exec A character string indicating the path to the GBLOCKS executable.
#' 
#' @details
#' Explanation of the routine taken from the Online Documentation: First, the
#' degree of conservation of every positions of the multiple alignment is
#' evaluated and classified as \emph{nonconserved}, \emph{conserved}, or
#' \emph{highly conserved}. All stretches of contiguous nonconserved positions
#' bigger than a certain value (\bold{b3}) are rejected. In such stretches,
#' alignments are normally ambiguous and, even when in some cases a unique
#' alignment could be given, multiple hidden substitutions make them inadequate
#' for phylogenetic analysis. In the remaining blocks, flanks are examined and
#' positions are removed until blocks are surrounded by highly conserved
#' positions at both flanks. This way, selected blocks are anchored by positions
#' that can be aligned with high confidence. Then, all gap positions -that can
#' be defined in three different ways (\bold{b5})- are removed. Furthermore,
#' nonconserved positions adjacent to a gap position are also eliminated until a
#' conserved position is reached, because regions adjacent to a gap are the most
#' difficult to align. Finally, small blocks (falling below the limit of
#' \bold{b4}) remaining after gap cleaning are also removed.
#' 
#' @returns A \code{matrix} of class \code{"DNAbin"}
#' @note 
#' From \code{phyloch} version 1.4-80 on, the defaults parameters of \code{gblocks} have been changed to be least conservative. \code{gblocks} was last updated and tested to work with Gblocks 0.91b. If you have problems getting the function to work with a newer version of Gblocks, please contact the package maintainer.
#' @references 
#' Castresana, J. 2000. Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. \emph{Molecular Biology and Evolution} \bold{17}, 540-552. 
#'   
#' Talavera, G., and J. Castresana. 2007. Improvement of phylogenies after removing divergent and ambiguously aligned blocks from protein sequence alignments. \emph{Systematic Biology} \bold{56}, 564-577. 
#'    
#' Gblocks website: \url{https://home.cc.umanitoba.ca/~psgendb/doc/Castresana/Gblocks_documentation.html}
#' @seealso
#' \code{\link{mafft}} and \code{\link{prank}} for sequence alignment; \code{\link{aliscore}} for another alignment masking algorithm.
#' @export

gblocks <- function(x, b1 = .5, b2 = b1, b3 = ncol(x), 
                    b4 = 2, b5 = "a", exec){
  
  if ( inherits(x, "alignment") ) stop("cannot handle class 'alignment'")
  if ( inherits(x, "list") ) stop("cannot handle unaligned sequences")
  
  ## check parameters:
  ## -----------------
  if ( b1 < .5 | b1 > 1 ) stop ("b1 not in [0.5, 1]")
  if ( b2 < b1 | b2 > 1 ) stop ("b2 not in [b1, 1]")
  if ( b3 < 0 | b4 > ncol(x) ) stop ("b3 not in [0, ", ncol(x), "]")
  if ( b4 < 2 | b4 > ncol(x) ) stop ("b4 not in [2, ", ncol(x), "]")
  b5 <- match.arg(b5, c("a", "h", "n"))
  
  if ( !file.exists(exec) ) 
    stop("executable '", exec, "' does not exist", sep = "")
  
  
  ntax <- nrow(x)
  b1 <- round(ntax * b1) + 1
  b2 <- round(ntax * b2) + 1
  
  cat("\n--- executing Gblocks ---")
  cat("\nminimum number of sequences for a conserved position :", b1)
  cat("\nminimum number of sequences for a flank position     :", b2)
  cat("\nmaximum number of contiguous nonconserved positions  :", b3)
  cat("\nminimum length of a block                            :", b4)
  cat("\nallowed gap positions                                :", b5)
  
  write.fas(x, "R2GBLOCK.fas")
  system(paste(exec, " R2GBLOCK.fas -t=d", 
               " -b1=", b1,
               " -b2=", b2,
               " -b3=", b3,
               " -b4=", b4,
               " -b5=", b5, 
               sep = "")
         )
  out <- read.fas("R2GBLOCK.fas-gb")
  unlink(list.files(pattern = "R2GBLOCK"))
  out
}