#' as.matrix.AAbin
#' @export
as.matrix.AAbin <- function(x, ...)
{
  if (is.matrix(x)) return(x)
  if (is.vector(x)) {
    dim(x) <- c(1, length(x))
    return(x)
  }
  s <- unique(lengths(x, use.names = FALSE))
  if (length(s) != 1)
    stop("AA sequences in list not of the same length.")
  n <- length(x)
  y <- matrix(raw(), n, s)
  for (i in seq_len(n)) y[i, ] <- x[[i]]
  rownames(y) <- names(x)
  class(y) <- "AAbin"
  y
}
