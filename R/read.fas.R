## This code is part of the ips package
## © C. Heibl 2014 (last update 2017-03-22)
#' @export

read.fas <- function(x, text){

  if (!missing(text)){
    x <- text
  } else {
    x <- scan(x, what = character(), quiet = TRUE)
  }

  if (!length(x)) stop("file is empty")

  start <- grep("^ {0,}>", x)

  ## parse taxon names
  ## -----------------
  h <- unlist(strsplit(x[start][1], ""))
  if (length(h) == 1){
    taxnames <- x[start + 1]
    space <- 2
  } else {
    taxnames <- x[start]
    taxnames <- gsub(">", "", taxnames)
    space <- 1
  }
  ntax <- length(taxnames)

  ## parse characters
  ## ----------------
  start <- c(start, length(x) + 1)
  obj <- vector("list", ntax)
  for (i in 1:ntax)
    obj[[i]] <- unlist(strsplit(gsub(" ", "", x[(start[i] + space):(start[i + 1] - 1)]), NULL))
  names(obj) <- taxnames

  unique_for_AA <- c("Q","E","I","L","F","P","U","O","J","Z","X","*")
  if(!any(unique(unlist(obj)) %in% unique_for_AA)){
    ## DNA sequences
    obj <- lapply(obj, tolower)
    obj <- as.DNAbin(obj)
  } else {
    ## AA sequences
    obj <- rpg::as.AAbin(obj)
  }
  if (length(unique(sapply(obj, length))) == 1)
    obj <- as.matrix(obj)
  return(obj)
}
