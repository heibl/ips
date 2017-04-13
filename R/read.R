#' @title Reading Sequence Files
#' @description Read DNA and amino acid sequences from FASTA, PHILIP, and NEXUS 
#'   formatted files.
#' @param x A character string, giving the file name.
#' @param text A character string in FASTA format.
#' @return An matrix (aligned sequences)  or list (unaligned sequences) of class
#'   \code{DNAbin} or \code{AAbin}.
#' @references Maddison, D.R., D.L. Swofford, and W.P. Maddison. 1997. NEXUS: an
#'   extensible file format for systematic information. \emph{Syst. Biol.} 
#'   \bold{46}: 590-621.
#' @seealso \code{\link{mafft}} and \code{\link{prank}} for sequence alignment,
#'   \code{\link{gblocks}} and \code{\link{aliscore}} for quality check and
#'   cleaning of sequence alignments, \code{\link{cbind.DNAbin}} for
#'   concatenation of sequence alignments.
#' @examples 
#' ## bark beetle COX1 sequences
#' data(ips.cox1)
#' ## create temporary file names
#' format <- c(".fas", ".phy", ".nex")
#' fn <- sapply(format, tempfile, 
#'              pattern = "ips", tmpdir = tempdir())
#' ## write sequences files
#' write.fas(ips.cox1, fn[".fas"])
#' write.phy(ips.cox1, fn[".phy"])
#' write.nex(ips.cox1, fn[".nex"])
#' ## read sequence files
#' fas <- read.fas(fn[".fas"])
#' phy <- read.phy(fn[".phy"])
#' nex <- read.nex(fn[".nex"])
#' ## remove sequence files
#'unlink(fn)
#' @name read

NULL