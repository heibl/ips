## This code is part of the ips package
## © C. Heibl 2014 (last update 2021-04-04)

#' @title Sequence Alignment with MAFFT
#' @description This function is a wrapper for MAFFT and can be used for 
#'   (profile) aligning of DNA and amino acid sequences.
#' @param x An object of class \code{DNAbin} or \code{AAbin}.
#' @param y An object of class \code{DNAbin} or \code{AAbin}, if given both
#'   \code{x} and \code{y} are preserved and aligned to each other ("profile
#'   alignment").
#' @param add A character string giving the method used for adding \code{y} to 
#'   \code{x}: \code{"add"}, \code{"addprofile"} (default), or any unambiguous 
#'   abbreviation of these.
#' @param method A character string giving the alignment method. Available 
#'   accuracy-oriented methods for less than 200 sequences are 
#'   \code{"localpair"}, \code{"globalpair"}, and \code{"genafpair"}; 
#'   \code{"retree 1"} and \code{"retree 2"} are for speed-oriented alignment. 
#'   The default is \code{"auto"}, which lets MAFFT choose an appropriate 
#'   alignment method.
#' @param maxiterate An integer giving the number of cycles of iterative 
#'   refinement to perform. Possible choices are \code{0}: progressive method, 
#'   no iterative refinement (default); \code{2}: two cycles of iterative 
#'   refinement; \code{1000}: at most 1000 cycles of iterative refinement.
#' @param op A numeric giving the \code{gap opening penalty} at group-to-group 
#'   alignment; default 1.53.
#' @param ep A numeric giving the offset value, which works like \code{gap 
#'   extension penalty}, for group-to-group alignment; default 0.0, but 0.123 is
#'   recommended if no long indels are expected.
#' @param gt An object of class \code{\link{phylo}} that is to be used as a
#'   guide tree during alignment. The default is to let MAFFT use its own
#'   internal guide tree
#' @param options A vector of mode character specifying additional arguments to 
#'   MAFFT, that are not included in \code{mafft} such as, e.g., 
#'   \code{--adjustdirection}.
#' @param thread Integer giving the number of physical cores MAFFT should use; 
#'   with \code{thread = -1} the number of cores is determined automatically.
#' @param exec A character string giving the path to the MAFFT executable 
#'   including its name, e.g. something like \code{/user/local/bin/mafft} under 
#'   UNIX-alikes.
#' @param quiet Logical, if set to \code{TRUE}, mafft progress is printed out on
#'   the screen.
#' @param file A character string indicating the filename of the output FASTA
#'   file; if this is missing the the alignment will be returned as  matrix of
#'   class \code{DNAbin} or \code{AAbin}.
#' @details \code{"localpair"} selects the \bold{L-INS-i} algorithm, probably 
#'   most accurate; recommended for <200 sequences; iterative refinement method 
#'   incorporating local pairwise alignment information.
#'   
#'   \code{"globalpair"} selects the \bold{G-INS-i} algorithm suitable for 
#'   sequences of similar lengths; recommended for <200 sequences; iterative 
#'   refinement method incorporating global pairwise alignment information.
#'   
#'   \code{"genafpair"} selects the \bold{E-INS-i} algorithm suitable for 
#'   sequences containing large unalignable regions; recommended for <200 
#'   sequences.
#'   
#'   \code{"retree 1"} selects the \bold{FFT-NS-1} algorithm, the simplest 
#'   progressive option in MAFFT; recommended for >200 sequences.
#'   
#'   \code{"retree 2"} selects the \bold{FFT-NS-2} algorithm that uses a second 
#'   iteration of alignment based on a guide tree computed from an FFT-NS-1 
#'   alignment; this is the default in MAFFT; recommended for >200 sequences.
#' @return A \code{matrix} of class \code{"DNAbin"} or \code{"AAbin"}.
#' @references Katoh, K. and H. Toh. 2008. Recent developments in the MAFFT 
#'   multiple sequence alignment program. \emph{Briefings in Bioinformatics} 
#'   \bold{9}: 286-298.
#'   
#'   Katoh, K., K.-i. Kuma, H. Toh, and T. Miyata. 2005. Mafft version 5: 
#'   improvement in accuracy of multiple sequence alignment. \emph{Nucleic Acids
#'   Research} \bold{33}: 511--518.
#'   
#'   Katoh, K., K. Misawa, K.-i. Kuma, and T. Miyata. 2002. Mafft: a novel 
#'   method for rapid multiple sequence alignment based on fast Fourier 
#'   transform. \emph{Nucleid Acids Research} \bold{30}: 3059--3066.
#'   
#'   \url{https://mafft.cbrc.jp/alignment/software/index.html}
#' @note \code{mafft} was last updated and tested to work with MAFFT 7.205. If 
#'   you have problems getting the function to work with a newer version of 
#'   MAFFT, please contact the package maintainer.
#' @seealso \code{\link{read.fas}} to import DNA sequences; \code{\link{prank}} 
#'   for another alignment algorithm; \code{\link{gblocks}} and 
#'   \code{\link{aliscore}} for alignment cleaning.
#' @importFrom phangorn write.phyDat
#' @export

mafft <- function(x, y, add, method = "auto", maxiterate = 0,
                  op = 1.53, ep = 0.0, gt = NULL, options,
                  thread = -1, exec, quiet, file){
  
  ## CHECKS and DEFINITIONS
  ## ----------------------
  if (missing(exec)) exec <- "/usr/local/bin/mafft"
  # info <- system2(command = exec, args = "--help", stdout = FALSE)

  if (!inherits(x, c("DNAbin", "AAbin")))
    stop("'x' is not of class 'DNAbin' or 'AAbin'")
  os <- .Platform$OS
  if (missing(quiet)) quiet <- TRUE
  qut <- ifelse(quiet, " --quiet ", " ")
  
  maxiterate <- match.arg(as.character(maxiterate), c("0", "2", "1000"))
  if (is.matrix(x)) x <- as.list(x)
  
  ## temporary input/output files
  ## ----------------------------
  fns <- vector(length = 4)
  for (i in seq_along(fns))
    fns[i] <- tempfile(pattern = "mafft", tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])
  method <- match.arg(method, c("auto", "localpair", "globalpair",
                                "genafpair", "parttree",
                                "retree 1", "retree 2"))
  
  ## guide tree
  ## ----------
  if (is.null(gt)){
    gt <- ""
  } else {
    if (!inherits(gt, "phylo"))
      stop("object \"gt\" is not of class \"phylo\"")
    if (!all(names(x) %in% gt$tip.label))
      stop("guide tree does not match sequence names")
    gt$tip.label <- match(names(x), gt$tip.label)
    if (is.null(gt$edge.length))
      gt$edge.length <- rep(1, nrow(gt$edge))
    phylo2mafft(gt, file = fns[4])
    gt <- paste(" --treein", fns[4], "")
  }
  
  ## multithreading
  ## --------------
  thread <- paste("--thread", thread)
  thread <- paste(rep(" ", 2), collapse = thread)
  
  ## additional arguments specified by the user
  ## ------------------------------------------
  if (missing(options)){
    options <- " "
  } else {
    options <- match.arg(options, c("--adjustdirection",
                                    "--adjustdirectionaccurately"))
    options <- paste(options, collapse = " ")
    options <- paste(rep(" ", 2), collapse = options)
  }
  
  ## write input files and prepare call to MAFFT
  ## -------------------------------------------
  if (missing(y)){
    if (inherits(x, "DNAbin")){ write.fas(x, fns[1])  }
    if (inherits(x, "AAbin")) { write.fas(x, fns[1]) }
    call.mafft <- paste(exec, " --", method, " --",
                        "maxiterate ", maxiterate, qut, "--op ", op,
                        " --ep ", ep, gt, options, thread,
                        fns[1], " > ", fns[3], sep = "")
  } else {
    if (!inherits(y, c("DNAbin", "AAbin"))) stop("'y' is not of class 'DNAbin' or 'AAbin'")
    if (is.matrix(y)) y <- as.list(y)
    if (missing(add)) add <- "addprofile"
    add <- match.arg(add, c("add", "addprofile"))
    add <- paste("--", add, sep = "")
    write.fas(x, fns[1])
    write.fas(y, fns[2])
    call.mafft <- paste(exec, qut, add, fns[2], fns[1], ">", fns[3])
  }
  if (!quiet) message(call.mafft)
  
  ## execute MAFFT on UNIX
  ## ---------------------
  if (os == "unix"){
    system(call.mafft, intern = FALSE, ignore.stdout = FALSE)
    res <- (file.info(fns[3])$size > 1) #allow line break in empty file
    if (res != 0) {
      res <- read.dna(fns[3],format="fasta")
    }
    ## execute MAFFT on WINDOWS
    ## ------------------------
  } else {
    res <- system(call.mafft, intern = TRUE, ignore.stderr = FALSE)
    if (length(grep("error|ERROR", res))){
      res <- 0
    }
    else {
      res <- read.dna(fns[3],format="fasta")
    }
  }
  unlink(fns[file.exists(fns)])
  if (!missing(file)){
    write.fas(res, file)
  } else {
    return(res)
  }
}
