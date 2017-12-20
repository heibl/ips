## This code is part of the ips package
## © C. Heibl 2014 (last update 2017-12-20)

#' @title Profile Alignment with MAFFT
#' @description Merge two or more DNA or amino acis sequence alignment by
#'   profile alignment with MAFFT.
#' @param subMSA A list of objects of class \code{"\link{DNAbin}"} or
#'   \code{"\link{AAbin}"}.
#' @param method A character string giving the alignment method. Available
#'   accuracy-oriented methods for less than 200 sequences are
#'   \code{"localpair"}, \code{"globalpair"}, and \code{"genafpair"};
#'   \code{"retree 1"} and \code{"retree 2"} are for speed-oriented alignment.
#'   The default is \code{"auto"}, which lets MAFFT choose an appropriate
#'   alignment method.
#' @param gt An object of class \code{\link{phylo}} that is to be used as a
#'   guide tree during alignment.
#' @param thread Integer giving the number of physical cores MAFFT should use;
#'   with \code{thread = -1} the number of cores is determined automatically.
#' @param exec A character string giving the path to the MAFFT executable
#'   including its name, e.g. something like \code{/user/local/bin/mafft} under
#'   UNIX-alikes.
#' @param quiet Logical, if set to \code{TRUE}, mafft progress is printed out on
#'   the screen.
#' @return An object of class \code{"\link{DNAbin}"} or \code{"\link{AAbin}"}.
#' @export

mafft.merge <- function(subMSA, method = "auto", gt,
                        thread = -1, exec, quiet = TRUE){

  quiet <- ifelse(quiet, "--quiet", "")

  ## set method
  ## ----------
  method <- match.arg(method, c("auto", "localpair", "globalpair",
                                "genafpair", "parttree",
                                "retree 1", "retree 2"))
  method <- paste("--", method, sep = "")

  ## guide tree
  ## ----------
  if (missing(gt)){
    gt <- ""
  } else {
    phylo2mafft(gt)
    gt <- "--treein tree.mafft"
  }

  ## create sub-MSA table
  ## --------------------
  n <- sapply(subMSA, nrow)
  subMSAtable <- vector(length = length(n))
  init <- 0
  for (i in seq_along(n)){
    nn <- 1:n[i] + init
    init <- max(nn)
    subMSAtable[i] <- paste(nn, collapse = " ")
  }

  ## prepare sequences input file
  ## ----------------------------
  subMSA <- lapply(subMSA, as.list)
  subMSA <- do.call(c, subMSA)
  names(subMSA) <- gsub("^.+[.]", "", names(subMSA))

  ## write input files
  ## -----------------
  fns <- vector(length = 3)
  for ( i in seq_along(fns) )
    fns[i] <- tempfile(pattern = "mafft",
                       tmpdir = tempdir()#,
                       #fileext = c(".txt", ".fas", ".fas")
                       )
  write(subMSAtable, fns[1])
  write.fas(subMSA, fns[2])

  ## assemble call to MAFFT
  ## ----------------------
  call.mafft <- paste(exec, method,
                      "--merge", fns[1],
                      quiet, gt, "--thread", thread,
                      fns[2], ">", fns[3])
#   cat(call.mafft)
#   if ( os == "unix" ){
    system(call.mafft, intern = FALSE, ignore.stdout = FALSE)
    res <- length(scan(fns[3], what = "c", quiet = TRUE))
    if (res != 0) res <- read.fas(fns[3])
#   }
  unlink(fns[file.exists(fns)])
  return(res)
}
