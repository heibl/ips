#' as.list.AAbin
#' @export
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



