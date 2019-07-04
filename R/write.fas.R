## This code is part of the ips package
## © C. Heibl 2014 (last update 2019-06-27)

#' @export

write.fas <- function(x, file, block.width = FALSE,
                      truncate = FALSE, append = FALSE){
  
  
  
  # x MUST be a list of class DNAbin
  # --------------------------------
  if (is.matrix(x)){
    x <- as.list(x)
  }
  if (is.data.frame(x)){
    y  <- as.list(x[, 1])
    names(y) <- rownames(x)
    x <- y
  }
  
  # taxonnames
  # ----------
  taxnames <- names(x)
  if (truncate){
    taxnames <- substring(taxnames, 1, truncate)
    if (any(duplicated(taxnames)))
      warning("truncation of taxon names created identical strings")
  }
  taxnames <- paste0(">", taxnames)
  
  # function to interleave sequences:
  # ---------------------------------
  create.block.widthd <- function(x, bp){
    x <- unlist(strsplit(x, ""))
    ncha <- length(x)
    if (!is.numeric(bp)) bp <- ncha
    nbpar <- ceiling(ncha / bp)
    out <- vector(length = nbpar)
    for (i in seq(along = out)){
      idmax <- i * bp
      if (idmax > ncha) idmax <- ncha
      idmin <- i * bp - bp + 1
      out[i] <- paste(x[idmin:idmax], collapse = "")
    }
    out
  }
  
  # assemble FASTA file:
  # --------------------
  s <- as.character(x)
  s <- lapply(s, paste, collapse = "")
  s <- lapply(s, create.block.widthd, bp = block.width)
  nbpar <- sapply(s, length)
  s <- unlist(s)
  fas <- vector(length = length(s) + length(taxnames))
  id <- 1
  for (i in seq(along = nbpar[-1]))
    id <- c(id, tail(id, 1) + nbpar[i] + 1)
  fas[id] <- taxnames
  fas[-id] <- unlist(s)
  
  # write FASTA file
  # ----------------
  if ( missing(file) ) {
    ## return character vector:
    return(fas)
  } else {
    if ( file == "" ) {
      ## print onto screen:
      cat(fas, sep = "\n")
    } else {
      ## write to file:
      write(fas, file = file, append = append)
    }
  }
}
