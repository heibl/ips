## This code is part of the ips package
## © C. Heibl 2015 (last update 2015-11-28)

formatSet <- function(x, arg = "lset"){
  
  x.p <- attr(x, "partition")
  x <- paste(names(x), unlist(x), sep = "=")
  if ( !is.null(x.p) ){
    x.p <- which(names(x) %in% x.p )
    x.p <- paste(x.p, collapse = ",")
    x.p <- paste("applyto=(", x.p, ")", sep = "")
    x <- c(x.p, x)
  }
  x <- paste(x, collapse = " ")
  x <- paste("\t", arg, " ", x, ";", sep = "")
  x
}