rc.core <- function(x){
  nn <- as.raw(c(136, 40, 72, 24, 240, 2, 4))
  names(nn) <- c("a", "c", "g", "t","n", "?", "-")
  a <- which(x %in% nn["a"])
  c <- which(x %in% nn["c"])
  g <- which(x %in% nn["g"])
  t <- which(x %in% nn["t"])
  x[a] <- nn["t"]
  x[c] <- nn["g"]
  x[g] <- nn["c"]
  x[t] <- nn["a"]
  rev(x)
}

rc <- function(x){
  
  for ( i in seq_along(x) )
    x[[i]] <- rc.core(x[[i]])
  x
}