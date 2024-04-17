## This code is part of the ips package
## © C. Heibl 2014 (last update 2019-11-27)

#' @rdname EmptyCells
#' @export

deleteEmptyCells <- function(DNAbin, margin = c(1, 2),
                             nset = c("-", "n", "?"),
                             quiet = FALSE){
  
  if ( !inherits(DNAbin, "DNAbin") ) 
    stop("'DNAbin' is not of class 'DNAbin'")
  
  ## convert character to raw
  
  ## IUPAC ambiguity code
  ## --------------------
  iupac <- c(n = 240, "?" = 2, "-" = 4,
    # a = 136, c = 40, g = 72, t = 24, 
    r = 192, y = 48, s = 96, w = 144, k = 80, m = 160, 
    b = 112, d = 208, h = 176, v = 224)
  nset <- iupac[nset]
  nset <- as.raw(nset)
  
  ## function that detects non-empty strings
  isNotEmpty <- function(x, nset){
    ifelse(all(unique(x) %in% nset), FALSE, TRUE)
  }
  
  size <- dim(DNAbin)
  
  ## rows  (margin == 1)
  if (1 %in% margin){
    rowind <- which(apply(DNAbin, 1, isNotEmpty, nset = nset))
    DNAbin <- DNAbin[rowind, ]
  }
  
  ## columns (margin == 2)
  if (2 %in% margin){
    colind <- which(apply(DNAbin, 2, isNotEmpty, nset = nset))
    DNAbin <- DNAbin[, colind]
  }
  
  ## screen output (if desired)
  if (!quiet) {
    size <- size - dim(DNAbin)
    rows <- ifelse(size[1] == 1, " row ", " rows ")
    cols <- ifelse(size[2] == 1, " column ", " columns ")
    message(size[1], rows, "deleted from alignment\n",
            size[2], cols, "deleted from alignment")
  }  
  DNAbin
}

