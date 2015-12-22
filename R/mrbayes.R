## This code is part of the ips package
## © C. Heibl 2014 (last update 2015-12-10)

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

