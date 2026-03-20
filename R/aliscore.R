## This code is part of the ips package
## © C. Heibl 2014 (last update 2019-06-27)


#' @title Masking of Sequence Alignments with ALISCORE
#' @description Provides a interface to \bold{Aliscore}, in order to remove
#'   problematic regions of a DNA sequence alignment.
#' @param x DNA sequences of class \code{DNAbin}.
#' @param gaps A vector of mode \code{"character"} indicating how gaps shall be
#'   treated: as \code{"5state"} or as \code{"ambiguous"}.
#' @param w An integer giving the size of the sliding window.
#' @param r An integer giving the number of random pairwise sequence
#'   comparisons; defaults to \code{4 * N}.
#' @param t \emph{Not yet implemented}.
#' @param l \emph{Not yet implemented}.
#' @param s \emph{Not yet implemented}.
#' @param o A vector of mode \code{"character"} containing outgroup taxon names.
#' @param exec A character string, giving the path to the Aliscore script.
#' @return A \code{matrix} of class \code{"DNAbin"}.
#' @note This function was developed with ALISCORE version 2.
#' @references Misof, B. and K. Misof. 2009. A Monte Carlo approach successfully identifies
#' randomness in multiple sequence alignments: a more objective means of data
#' exclusion. \emph{Syst. Biol.} \bold{58}: 21--34.
#' @references Kueck, P., K. Meusemann, J. Dambach, B. Thormann, B.M. von Reumont, J.W.
#' Waegele and B. Misof. 2010. Parametric and non-parametric masking of
#' randomness in sequence alignments can be improved and leads to better
#' resolved trees. \emph{Frontiers in Zoology} \bold{7}: 10.
#' @seealso \code{\link{mafft}} and \code{\link{prank}} for multiple sequence
#'   alignment; \code{\link{gblocks}} for another alignment masking algorithm.
#' @examples
#' data(ips.28S)
#' \dontrun{aliscore(ips.28S)}
#' @export

aliscore <- function(x, gaps = "5state", w = 6, r, t, l, s, o, exec){
	
  ## Check path to executable
  ## ------------------------
  if (missing(exec)){
    ## Try to guess location of executable
    exec <- list.files(path = "/Applications", pattern = "(Aliscore)|(ALISCORE)")
    exec <- list.files(path = file.path("/Applications", exec), recursive = TRUE,
                       pattern = "Aliscore.02.2.pl", full.names = TRUE)
    if (!length(exec)) stop("path to executable not given")
    message("Using '", exec, "'", appendLF = TRUE)
  } else {
    if (!file.exists(exec)) 
      stop("executable '", exec, "' does not exist", sep = "")
  }
  
  in_file <- tempfile("input.fas")
	write.fas(x, in_file)
  
  ## Prepare call of Perl script
  ## ---------------------------
  N <- ifelse(gaps == "5state", "", "-N") # treatment of gaps
  w <- paste("-w", w) # window size
  r <- ifelse(missing(r), "", paste("-r", r ))
  if (!missing(t)) stop("option -t not yet implemented")
	if (!missing(l)) stop("option -l not yet implemented")
	if (!missing(s)) stop("option -s not yet implemented")
	o <- ifelse(missing(o), "", paste("-o", paste(o, collapse = ",") ))
	call <- paste("perl", exec, "-i", in_file, N, w, r, o)
	
	## Need to change @INC to locate Aliscore_module
	## ---------------------------------------------
	atINC <- basename(exec)
	atINC <- gsub(atINC, "", exec)
	atINC <- paste0("export PERL5LIB=", atINC, ";")
	
	## Execute, parse results and delete temporary files
	## ----------------------------------------
	system(paste(atINC, call))
	id <- scan(paste0(in_file, "_List_l_all.txt"),
             sep = " ", quiet = TRUE)
	# unlink(tempdir(), recursive = TRUE)
	id <- as.numeric(id)
	
	## Delete unreliable columns from MSA
	## ----------------------------------
	nid <- length(id)
	if (!nid) {
	  message("\nALISCORE did not remove any characters")
	} else {
	  x <- x[, -id]
	  message("\nALISCORE removed ", nid, " characters")
	}
	x
}
