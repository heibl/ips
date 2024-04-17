## This code is part of the ips package
## © C. Heibl 2014 (last update 2018-01-22)

#' @rdname read
#' @export

read.fas <- function(x, text){

  if (!missing(text)){
    x <- text
  } else {
    x <- scan(x, what = character(), quiet = TRUE)
  }

  if (!length(x)) stop("file is empty")

  start <- grep("^ {0,}>", x)
  
  ## Parse taxon names
  ## -----------------
  h <- unlist(strsplit(x[start][1], ""))
  if (length(h) == 1){
    taxnames <- x[start + 1]
    nlines_metadata <- 2
  } else {
    ## Metadata that contain spaces will be broken into several lines
    ## and we will have to concatenata them (e.g. Fasta from ZSM)
    ## ----------------------------------------------------------
    nlines_metadata <- diff(start[1:2]) - 1
    if (nlines_metadata > 2){
      id <- lapply(start, function(start, n) start + (1:n) - 1, n = nlines_metadata)
      taxnames <- sapply(id, function(id, x) paste(x[id], collapse = "_"), x = x)
    } else {
      taxnames <- x[start]
    }
    taxnames <- gsub(">", "", taxnames)
  }
  ntax <- length(taxnames)

  ## Parse characters
  ## ----------------
  start <- c(start, length(x) + 1)
  obj <- vector("list", ntax)
  for (i in 1:ntax){
    obj[[i]] <- unlist(strsplit(gsub(" ", "", x[(start[i] + nlines_metadata):(start[i + 1] - 1)]), NULL))
  }
  names(obj) <- taxnames

  ## determine if sequences are DNA or AA:
  ## these strings are unique to amino acid sequences
  aa_string <- c("Q","E","I","L","F","P","U","O","J","Z","X","*")
  if (any(toupper(unlist(obj))  %in%  aa_string)){
    ## AA sequences
    obj <- as.AAbin(obj)
  } else {
    ## DNA sequences
    obj <- lapply(obj, tolower)
    obj <- as.DNAbin(obj)
  }
  if (length(unique(sapply(obj, length))) == 1)
    obj <- as.matrix(obj)
  return(obj)
}
