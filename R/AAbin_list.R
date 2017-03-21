## Franz Krah [2017-03-21]
## These are simple add-on functions to R package ape
## Ape currently does not support lists and bit-level coding scheme
## for amino acid sequences
## The functions here implement the former, not the latter.

as.list.AAbin <- function(x, ...)
{
  if (is.list(x)) return(x)
  if (is.null(dim(x))) obj <- list(x) # cause is.vector() doesn't work
  else { # matrix
    n <- nrow(x)
    obj <- vector("list", n)
    for (i in seq_len(n)) obj[[i]] <- x[i, , drop = TRUE]
    names(obj) <- rownames(x)
  }
  class(obj) <- "AAbin"
  obj
}


c.AAbin <- function(..., recursive = FALSE)
{
  if (!all(unlist(lapply(list(...), is.list))))
    stop("the 'c' method for \"AAbin\" accepts only lists")
  structure(NextMethod("c"), class = "AAbin")
}


as.AAbin.list <- function(x, ...)
{
  obj <- lapply(x, as.AAbin)
  class(obj) <- "AAbin"
  obj
}
