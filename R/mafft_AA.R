## This code is part of the ips package
## © C. Heibl 2014 (last update 2017-03-15)

#' @export

mafft_AA <- function(x, y, add, method = "auto", maxiterate = 0,
                  op = 1.53, ep = 0.0, gt, options,
                  thread = -1, exec, quiet){

  ## CHECKS and DEFINITIONS
  ## ----------------------
  if (!inherits(x, c("DNAbin", "AAbin")))
    stop("'x' is not of class 'DNAbin' or 'AAbin'")
  os <- .Platform$OS
  if (missing(quiet)) quiet <- TRUE
  qut <- ifelse(quiet, " --quiet ", " ")
  if (missing(exec)) exec <- "/usr/local/bin/mafft"
  maxiterate <- match.arg(as.character(maxiterate), c("0", "2", "1000"))

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
  if (missing(gt)){
    gt <- ""
  } else {
    if (!inherits(gt, "phylo"))
      stop("object \"gt\" is not of class \"phylo\"")
    if (!all(names(x) %in% gt$tip.label))
      stop("guide tree does not match sequence names")
    gt$tip.label <- match(names(x), gt$tip.label)
    if (!is.binary.tree(gt))
      gt <- multi2di(gt)
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
    if (inherits(x, "AAbin")) { write.phyDat(x, file = fns[1], format="fasta") }
    call.mafft <- paste(exec, " --", method, " --",
                        "maxiterate ", maxiterate, qut, "--op ", op,
                        " --ep ", ep, gt, options, thread,
                        fns[1], " > ", fns[3], sep = "")
  } else {
	  if (!inherits(y, "DNAbin")) stop("'y' is not of class 'DNAbin'")
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
    res <- length(scan(fns[3], what = "c", quiet = TRUE))
    if (res != 0) {
      #res <- read.fas(fns[3])
      if (inherits(x, "DNAbin")) { res <- read.fas(fns[3]) }
      if (inherits(x, "AAbin" )) { res <- seqinr::read.fasta(fns[3], seqtype  ="AA")   }
    }

  ## execute MAFFT on WINDOWS
  ## ------------------------
  } else {
    res <- system(call.mafft, intern = TRUE, ignore.stderr = FALSE)
    if (length(grep("error|ERROR", res))){
      res <- 0
    }
    else {
      # res <- read.fas(fns[3])
      if (inherits(x, "DNAbin")) { res <- read.fas(fns[3]) }
      if (inherits(x, "AAbin" )) { res <- seqinr::read.fasta(fns[3], seqtype  ="AA")}
    }
  }
  unlink(fns[file.exists(fns)])
  res
}

