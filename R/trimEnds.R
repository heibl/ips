## This code is part of the ips package
## © C. Heibl 2014 (last update 2019-11-05)

#' @title Trim Alignment Ends
#' @description Trims both ends of a DNA sequence alignment to the first and
#'   last alignment positions that contain a minimum number of IUPAC base
#'   characters (\code{"a"}, \code{"c"}, \code{"g"}, \code{"t"}, \code{"r"},
#'   \code{"y"}, \code{"s"}, \code{"w"}, \code{"k"}, \code{"m"}, \code{"b"},
#'   \code{"d"}, \code{"h"}, \code{"v"}). In addition, all gap characters
#'   (\code{"-"}) beyond the first and last base characters of each sequence are
#'   replaced by the  character \code{"n"}.
#' @param x An object of class \code{DNAbin}.
#' @param min.n.seq A \code{numeric} giving the required minimum number of
#'   sequences having an non-ambiguous base character (a, c, g, t) in the first
#'   and last position of the alignment; defaults to \code{4}, which is the
#'   minimum number of sequences needed to produce a non-trivial unrooted
#'   topology. Can also be given as a fraction.
#' @return An object of class \code{DNAbin}.
#' @seealso \code{\link{deleteEmptyCells}}, \code{\link{deleteGaps}}
#' @examples
#' # simple example alignment:
#' x <- structure(list(nb = 5, seq = c("acaaggtaca", "-caaggtac-",
#' "acaaggtaca", "aca--gtaca", "-ccaggta--"), nam = LETTERS[1:5]),
#' .Names = c("nb", "seq", "nam"), class = "alignment")
#' # convert to DNAbin:
#' x <- as.DNAbin(x)
#' # fill missing nucleotides:
#' x <- trimEnds(x)
#' # show results:
#' as.character(x[2, ])
#' @export

trimEnds <- function(x, min.n.seq = 4){
  
  if ( !inherits(x, "DNAbin") ){
    stop("'x' is not of class 'DNAbin'")
  }
  if ( !is.matrix(x) ){
    stop("'x' must be a matrix")
  }
  
  ## Store confidence stores; if not present cs == NULL
  ## --------------------------------------------------
  cs <- attr(x, "cs")
  
  ## Turn fraction into numbers
  ## --------------------------
  if (min.n.seq < 1){
    min.n.seq <- ceiling(nrow(x) * min.n.seq)
  }
  
  ## If alignment has less then min.n.seq sequences,
  ## min.n.seq has to be adjusted
  ## ----------------------------
  min.n.seq <- min(nrow(x), min.n.seq)
  
  ## Replace terminal '-' with 'N'
  ## -----------------------------
  replaceWithN <- function(x){
    
    n <- vector()
    
    ## head (5'-end)
    id <- which(x == as.raw(4))
    if ( 1 %in% id ) n <- c(n, which(id == 1:length(id)))
    
    ## tail (3'-end)
    id <- which(rev(x == as.raw(4)))
    if ( 1 %in% id ) n <- c(n, (length(x):1)[which(id == 1:length(id))])
    
    ## replace - by N
    if ( length(n) > 0 ){
      x[n] <- as.raw(240)
    }
    x
  }
  x <- t(apply(x, 1, replaceWithN))
  class(x) <- "DNAbin"
  
  ## Remove 'sandspit' pattern

  ## -------------------------
  removeSandspit <- function(x){
    
    ## anything to do?
    id <- which(x == as.raw(4))
    if (length(id) == 0) return(x)
    
    ## head (5'-end)
    n <- vector()
    lagoon <- id[id == min(id) + (1:length(id)) - 1]
    spit <- 1:(min(lagoon) - 1)
    if ( length(spit) <= 10 & length(lagoon) >= 5 ){
      n <- c(n, union(spit, lagoon))
    }
    
    ## tail (3'-end)
    id <- which(rev(x == as.raw(4)))
    lagoon <- lagoon[lagoon == min(lagoon) + (1:length(lagoon)) - 1]
    spit <- 1:(min(lagoon) - 1)
    if ( length(spit) <= 10 & length(lagoon) >= 5 ){
      n <- c(n, (length(id):1)[union(spit, lagoon)])
    }
    if ( length(n) > 0 ){
      x[n] <- as.raw(240)
    }
    x
  }
  x <- t(apply(x, 1, removeSandspit))
  class(x) <- "DNAbin"
  
  ## Trim ends to 'min.n.seq' bases
  ## ------------------------------
  iupac <- c(a = 136, c = 40, g = 72, t = 24, 
             r = 192, y = 48, s = 96, w = 144, k = 80, m = 160, 
             b = 112, d = 208, h = 176, v = 224)
  iupac <- as.raw(iupac)
  percentInformation <- function(x, iupac){
    length(which(x %in% iupac))
  }
  m <- apply(x, 2, percentInformation, iupac)
  if ( max(m) < min.n.seq ) stop("alignment contains less sequences then required")
  m <- range(which(m >= min.n.seq))
  m <- seq(from = m[1], to = m[2])
  x <- x[, m]
  
  ## Trim and reappend confidence scores
  ## -----------------------------------
  if (!is.null(cs)){
    if (is.matrix(cs)){
      cs <- cs[, m]
    } else {
      cs <- cs[m]
    }
    attr(x, "cs") <- cs
  }
  x
}



