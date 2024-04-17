## This code is part of the ips package
## © C. Heibl 2014 (last update 2017-03-14)

#' @export

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
