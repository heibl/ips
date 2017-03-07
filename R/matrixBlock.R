## This code is part of the ips package
## © C. Heibl 2014 (last update 2016-01-26)

matrixBlock <- function(x, block.width){
  
  # indices of partitions
  # ---------------------
  pt <- seq(from = 1, to = ncol(x), by = block.width)
  pt <- data.frame(from = pt, to = c(pt[-1] - 1, ncol(x)))
  
  # assemble matrix
  # ---------------
  if ( inherits(x, "DNAbin") ){
    m <- as.character(x)
    m <- apply(m, 1, paste, collapse = "")
    m <- apply(pt, 1, function(m, i) substr(m, i[1], i[2]), m = m)
    m <- as.vector(m)
    m <- paste(rep(rownames(x), nrow(pt)), unlist(m))
  } else {
    m <- apply(pt, 1, function(m, i) m[, i[1]:i[2]], m = x)
    m <- lapply(m, function(z) apply(z, 1, paste, collapse = ""))
    m <- unlist(m)
    m <- paste(rep(rownames(x), nrow(pt)), unlist(m))
  }
  m
}