## This code is part of the ips package
## © C. Heibl 2014 (last update 2016-08-01)

trimEnds <- function(x, min.n.seq = 4){
  
  if ( !inherits(x, "DNAbin") ){
    stop("'x' is not of class 'DNAbin'")
  }
  if ( !is.matrix(x) ){
    stop("'x' must be a matrix")
  }
  
  ## replace terminal '-' with 'N'
  ## -----------------------------
  replaceWithN <- function(x){
    
    n <- vector()
    
    ## head (5'-end)
    id <- which(x == as.raw(4))
    if ( 1 %in% id ) n <- c(n, which(id == 1:length(id)))
    
    ## tail (3'-end)
    id <- which(rev(x == as.raw(4)))
    if ( 1 %in% id ) n <- c(n, (length(x):1)[which(id == 1:length(id))])
    
    ## replace - by N
    if ( length(n) > 0 ){
      x[n] <- as.raw(240)
    }
    x
  }
  x <- t(apply(x, 1, replaceWithN))
  class(x) <- "DNAbin"
  
  ## remove 'sandspit' pattern
  ## -------------------------
  removeSandspit <- function(x){
    
    ## anything to do?
    id <- which(x == as.raw(4))
    if ( length(id) == 0 ) return(x)
    
    ## head (5'-end)
    n <- vector()
    lagoon <- id[id == min(id) + (1:length(id)) - 1]
    spit <- 1:(min(lagoon) - 1)
    if ( length(spit) <= 10 & length(lagoon) >= 5 ){
      n <- c(n, union(spit, lagoon))
    }
    
    ## tail (3'-end)
    id <- which(rev(x == as.raw(4)))
    lagoon <- lagoon[lagoon == min(lagoon) + (1:length(lagoon)) - 1]
    spit <- 1:(min(lagoon) - 1)
    if ( length(spit) <= 10 & length(lagoon) >= 5 ){
      n <- c(n, (length(id):1)[union(spit, lagoon)])
    }
    if ( length(n) > 0 ){
      x[n] <- as.raw(240)
    }
    x
  }
  x <- t(apply(x, 1, removeSandspit))
  class(x) <- "DNAbin"
  
  ## trim ends to 'min.n.seq' bases
  ## ------------------------------
  b <- as.raw(c(136, 40, 72, 24))
  percentInformation <- function(x, b){
    length(which(x %in% b))
  }
  m <- apply(x, 2, percentInformation, b)
  if ( max(m) < min.n.seq ) stop("alignment contains less sequences then required")
  m <- range(which(m >= min.n.seq))
  m <- seq(from = m[1], to = m[2])
  x <- x[, m]
  x
}



