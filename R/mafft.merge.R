## This code is part of the ips package
## © C. Heibl 2014 (last update 2017-04-12)

#' @export

mafft.merge <- function(subMSA, method = "auto", gt,
                        quiet = TRUE, thread = -1, exec){

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
