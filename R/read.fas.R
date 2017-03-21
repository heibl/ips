## This code is part of the ips package
## © C. Heibl 2014 (last update 2017-03-17)

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

  dna_string <- c("n", "?", "-", "r", "y", "s", "w", "k", "m", "b", "d", "h", "v")
  if (!all(tolower(unlist(obj))  %in%  dna_string)){
    ## DNA sequences
    obj <- lapply(obj, tolower)
    obj <- as.DNAbin(obj)
  } else {
    ## AA sequences
    obj <- as.AAbin(obj)

    ## CH [2017-03-17]
    ## Hier gehts im Moment nicht weiter, da AAbin keine Listen unterstützt.
    ## phangorn::read.aa gibt eine Fehlermeldung:
    ## Fehler in phangorn::read.aa(x) :
    ##   the first line of the file must contain the dimensions of the data
    ## Es wäre wohl sinnvoll die AAbin Klasse zu erweitern

    ## FK [2017-03-21]
    ## habe die AAbin Klasse erweitert. Liegt aktuell im File 'AAbin_list.R'
    ## Sollen wir das an Emmanuel schicken?
  }
  if (length(unique(sapply(obj, length))) == 1)
    obj <- as.matrix(obj)
  return(obj)
}
