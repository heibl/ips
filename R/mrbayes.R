## This code is part of the ips package
## © C. Heibl 2014 (last update 2025-08-30)

#' @title Bayesian MCMC Tree Search with MrBayes
#' @description A wrapper for Bayesian phylogenetic tree search through MrBayes
#'   (Ronquist & Huelsenbeck, 2003) with either DNA (\code{mrbayes}) or
#'   morphological (\code{mrbayes.mixed}) data.
#' @param x An object of class \code{\link[ape]{DNAbin}} in the case of
#'   \code{mrbayes} or a matrix of mode character in the case of
#'   \code{mrbayes.mixed}.
#' @param file A character string, giving the name of the MrBayes input file.
#' @param lset The output of a call to \code{\link{mrbayes.lset}}.
#' @param prset The output of a call to \code{\link{mrbayes.prset}}.
#' @param mcmc The output of a call to \code{\link{mrbayes.mcmc}}.
#' @param unlink xxx
#' @param constraint xxx
#' @param burnin An integer, the number of samples from the MCMC to be discarded
#'   prior to further analysis.
#' @param contype A character string; the type of consensus tree calculated from
#'   the posterior distribution of trees: either \code{"halfcompat"}
#'   (majority-rule consensus tree) or \code{"allcombat"} (strict consensus
#'   tree).
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

mrbayes <- function(x, file = "", lset, prset, mcmc, unlink, constraint,
                    burnin = 10, contype = "allcompat", run = FALSE){
  
  if ( !is.list(x) ) x <- list(x)
  
  ## lset:
  if ( !missing(lset) ){
    if ( !all(sapply(lset, is.list)) ){
      lset <- list(lset)
    }
    lset <- sapply(lset, formatSet, arg = "lset")
  } else {
    lset <- NULL
  }
  
  ## prset:
  if ( !missing(prset) ){
    if ( !all(sapply(prset, is.list)) ){
      prset <- list(prset)
    }
    prset <- sapply(prset, formatSet, arg = "prset")
  } else {
    prset <- NULL
  }
 
  ## mcmc:
  if ( !missing(mcmc) ){
    mcmc <- paste(names(mcmc), unlist(mcmc), sep = "=")
    mcmc <- paste(mcmc, collapse = " ")
    mcmc <- paste("\tmcmc ", mcmc, ";", sep = "")
  } else {
    mcmc <- "\tmcmc;"
  }

  ## handle partitions
  ## -----------------
  if ( length(x) > 1 ){
    
    ## get position of nucleotides:
    p <- cbind(rep(1, length(x)), sapply(x, ncol))
    if ( nrow(p) > 1 ) {
      for ( i in 2:nrow(p) ){
        p[i, 1] <- p[i - 1, 2] + 1
        p[i, 2] <- p[i, 1] + p[i, 2] -1
      }
    }
    partition <- c(
      paste("\tcharset ", rownames(p), " = ", 
          apply(p, 1, paste, collapse = "-"), ";", 
          sep = ""),
      paste("\tpartition a = ", length(x), ": ", 
            paste(rownames(p), collapse = ", "),
            ";", sep = ""),
      "\tset partition = a;"
    ) 
  } else {
    partition <- NULL
  }
  
  ## constraints
  ## -----------
  if ( !missing(unlink) ){
    unlink <- lapply(unlink, function(z) paste("(", z, ")", sep = ""))
    unlink <- paste(names(unlink), "=", unlink, collapse = " ")
    unlink <- paste("\tunlink ", unlink, ";", sep = "")
  } else {
    unlink <- NULL
  }
  
  ## constraints
  ## -----------
  if ( !missing(constraint) ){
    constraint <- lapply(constraint, paste, collapse = " ")
    constraint <- paste("\tconstraint ", names(constraint), " = ", 
                        constraint, ";", sep = "")
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
    paste("\tsumt filename=", file, " burnin=", 		
          burnin, " contype=", contype, ";", sep = ""),
    "end;"
  )
  
  # assemble and print NEXUS file
  # -----------------------------
  nexus <- write.nex(x, block.width = FALSE)
  nexus <- c(nexus, bayes)
  if ( file == "" ){
    cat(nexus, sep = "\n")
  } else {
    write(nexus, file = file)
  }
  
  # start mrbayes
  # -------------
  if ( run ){
    if ( .Platform$OS.type == "unix" ){
      system(paste("mb > execute", file))
    } else {
      system(paste("mrbayes", file)) # bug fix by Klaus Schliep (2015-10-27)
    }												
    tr <- read.nexus(paste(file, ".con.tre", sep = ""))
    tr
  }
}

