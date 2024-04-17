## This code is part of the ips package
## © C. Heibl 2016 (last update 2016-11-23)

#' @title Unlist To First Level Only
#' @description Does the same as \code{unlist}, but recurses only
#' one level.
#' @param z A list of lists.
#' @param use.names Logical, indicating if element names from the
#' element should be preserved.
#' @export

unlistFirstLevel <- function(z, use.names = TRUE){
  
  ## get ID of elements which are list ...
  id <- sapply(z, is.list)
  ## ... but not of class 'phylo'
  id <- id & !sapply(z, inherits, what =  "phylo")
  
  if ( any(id) ){
    
    ## unlist
    zz <- vector(mode = "list")
    for ( i in seq_along(z) ){
      if ( id[i] ) zz <- c(zz, z[[i]])
      else zz <- c(zz, z[i])
    }
    
    ## preserve names
    if ( use.names ){
      d <- sapply(z, length)
      names(zz) <- rep(names(d), d)
    }
    
    z <- zz
  }
  z
}