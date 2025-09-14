## This code is part of the ips package
## © C. Heibl 2014 (last update 2025-08-30)

#' @title Bayesian MCMC Tree Search with MrBayes
#' @description Provides a wrapper for Bayesian phylogenetic tree search through
#'   MrBayes (Huelsenbeck & Ronquist, 2001; Ronquist & Huelsenbeck, 2003).
#' @param x An object of class \code{\link[ape]{DNAbin}} in the case of
#'   \code{mrbayes} or a matrix of mode character in the case of
#'   \code{mrbayes.mixed}.
#' @param file A character string, giving the name of the MrBayes input file.
#' @param lset A list as returned by \code{\link{mrbayes.prset}} containing the
#'   parameter setting for the prior distributions.
#' @param prset A list as returned by \code{\link{mrbayes.prset}} containing the
#'   parameter setting for the prior distributions.
#' @param mcmc A list as returned by \code{\link{mrbayes.mcmc}} containing the
#'   parameter setting for the Markov chain Monte Carlo (MCMC).
#' @param unlink xxx
#' @param constraint xxx
#' @param burnin An integer; the number of samples from the MCMC to be discarded
#' @param exec A character string giving the full path of the MrBayes program.
#' @param run Logical; \code{run = FALSE} will only print the NEXUS file,
#'   \code{run = TRUE} will also start the MCMC runs, if \code{exec} is
#'   correctly specified.
#' @details \code{mrbayes} was last updated and tested with MrBayes
#'   \bold{v3.2.2} under R 3.1.0 on a x86_64-apple-darwin10.8.0 (64-bit)
#'   platform. It is intended to offer a simply parameterized building block for
#'   larger scripts.
#' @return None; a NEXUS file with MrBayes block is written to a file and, if
#'   \code{run = TRUE}, the MCMC runs in MrBayes are started.
#' @seealso \code{\link{mafft}} and \code{\link{prank}} for sequence alignment;
#'   \code{\link{raxml}} for maximum likelihood tree search.
#' @references 
#' J. P. Huelsenbeck & Ronquist F. 2001. MrBayes: Bayesian inference of phylogenetic trees. \emph{Bioinformatics} \bold{17}: 754-755.

#' Ronquist F. & J. P. Huelsenbeck. 2003. MrBayes 3: Bayesian phylogenetic inference under mixed models. \emph{Biometrics} \bold{19}: 1572-1574.

#' MrBayes website: \url{https://mrbayes.sourceforge.net/}.
#' @examples
#' data(ips.cox1)
#' x <- ips.cox1[, 100:140] # tiny alignment
#' mrbayes(x, file = "", mcmc = mrbayes.mcmc(ngen = 100), run = FALSE)
#' \dontrun{
#' library(phangorn)
#' tree <- rtree(10)
#' Y1 <- simSeq(tree, l = 20)
#' Y2 <- simSeq(tree, l = 20, type = "USER", levels=c("0", "1"))
#' Y <- cbind(as.character(Y1), as.character(Y2))
# 'mrbayes.mixed(Y, file = "", ngen = 100, run = FALSE)
#' }
#' @export 




=======
#' @param run Logical, \code{run = FALSE} will only print the NEXUS file,
#'   \code{run = TRUE} will also start the MCMC runs, if the \code{path}
#'   argument is correctly specified.
#' @details
#' \code{mrbayes} was last updated and tested with MrBayes \bold{v3.2.2} under R
#' 3.1.0 on a x86_64-apple-darwin10.8.0 (64-bit) platform. It is intended to
#' offer a simply parameterized building block for larger scripts.
#' @returns None, a NEXUS file with MrBayes block is written to a file and, if
#'   \code{run = TRUE}, the MCMC runs in MrBayes are started.
#' @references 
#' J. P. Huelsenbeck & Ronquist F. 2001. MrBayes: Bayesian inference of phylogenetic trees. \emph{Bioinformatics} \bold{17}: 754-755.
#'
#' Ronquist F. & J. P. Huelsenbeck. 2003. MrBayes 3: Bayesian phylogenetic inference under mixed models. \emph{Biometrics} \bold{19}: 1572-1574.
#' 
#' MrBayes website: \url{https://nbisweden.github.io/MrBayes/}.
#' @seealso \code{\link{mafft}} and \code{\link{prank}} for sequence alignment; \code{\link{raxml}} for maximum likelihood tree search.
#' @examples
#' data(ips.cox1)
#' x <- ips.cox1[, 100:140] # tiny alignment
#' # mrbayes(x, file = "", ngen = 100, run = FALSE)
#' 
#' \dontrun{
#'   
#'   library(phangorn)
#'   tree <- rtree(10)
#'   Y1 <- simSeq(tree, l = 20)
#'   Y2 <- simSeq(tree, l = 20, type = "USER", levels=c("0", "1"))
#'   Y <- cbind(as.character(Y1), as.character(Y2))
#'   mrbayes.mixed(Y, file = "", ngen = 100, run = FALSE)
#' }
#' 
#' @importFrom ape read.nexus
#' @export
>>>>>>> kurt

mrbayes <- function(x, file = "", lset, prset, mcmc, unlink, constraint,
                    burnin = 10, contype = "allcompat", exec, run = FALSE){

  if (!is.list(x)) x <- list(x)

  ## lset:
  if (!missing(lset)){
    if (!all(sapply(lset, is.list))){
      lset <- list(lset)
    }
    lset <- sapply(lset, formatSet, arg = "lset")
  } else {
    lset <- NULL
  }

  ## prset:
  if (!missing(prset)){
    if (!all(sapply(prset, is.list))){
      prset <- list(prset)
    }
    prset <- sapply(prset, formatSet, arg = "prset")
  } else {
    prset <- NULL
  }

  ## mcmc:
  if (!missing(mcmc)){
    mcmc <- paste(names(mcmc), unlist(mcmc), sep = "=")
    mcmc <- paste(mcmc, collapse = " ")
    mcmc <- paste0("\tmcmc ", mcmc, ";")
  } else {
    mcmc <- "\tmcmc;"
  }

  ## handle partitions
  ## -----------------
  if (length(x) > 1){

    ## get position of nucleotides:
    p <- cbind(rep(1, length(x)), sapply(x, ncol))
    if (nrow(p) > 1) {
      for ( i in 2:nrow(p) ){
        p[i, 1] <- p[i - 1, 2] + 1
        p[i, 2] <- p[i, 1] + p[i, 2] -1
      }
    }
    partition <- c(
      paste0("\tcharset ", rownames(p), " = ",
          apply(p, 1, paste, collapse = "-"), ";"),
      paste0("\tpartition a = ", length(x), ": ",
            paste(rownames(p), collapse = ", "), ";"),
      "\tset partition = a;"
    )
  } else {
    partition <- NULL
  }

  ## constraints
  ## -----------
  if (!missing(unlink)){
    unlink <- lapply(unlink, function(z) paste0("(", z, ")"))
    unlink <- paste(names(unlink), "=", unlink, collapse = " ")
    unlink <- paste0("\tunlink ", unlink, ";")
  } else {
    unlink <- NULL
  }

  ## constraints
  ## -----------
  if (!missing(constraint)){
    constraint <- lapply(constraint, paste, collapse = " ")
    constraint <- paste0("\tconstraint ", names(constraint), " = ",
                        constraint, ";")
  } else {
    constraint <- NULL
  }

  # create MrBayes block
  # --------------------
  bayes <- c(
    "\nbegin mrbayes;",
    partition,
    constraint,
    lset,
    unlink,
    prset,
    mcmc,
    paste0("\tsumt filename=", file, " burnin=",
          burnin, " contype=", contype, ";"),
    "end;"
  )

  # assemble and print NEXUS file
  # -----------------------------
  nexus <- write.nex(x, block.width = FALSE)
  nexus <- c(nexus, bayes)
  if (file == ""){
    cat(nexus, sep = "\n")
  } else {
    write(nexus, file = file)
  }

  # start mrbayes
  # -------------
  if (run){

    ## localize executable
    ## -------------------
    if (missing(exec)){
      exec <- ifelse(.Platform$OS.type == "unix", "mb", "mrbayes")
    }

    ## execute MrBayes
    ## ---------------
    if (.Platform$OS.type == "unix"){
      mb_call <- paste(exec, "> execute", file)
    } else {
      mb_call <- paste(exec, file) # bug fix by Klaus Schliep (2015-10-27)
    }
    system(mb_call)
    tr <- read.nexus(paste0(file, ".con.tre"))
    tr
  }
}
