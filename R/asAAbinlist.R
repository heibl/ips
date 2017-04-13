#' as.AAbin.list
#' @export

as.AAbin.list <- function(x, ...)
{
  obj <- lapply(x, as.AAbin)
  class(obj) <- "AAbin"
  obj
}
