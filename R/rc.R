## This code is part of the ips package
## © C. Heibl 2018 (last update 2022-02-05)

#' @title Reverse-Complement of DNA sequences
#' @description Reverse, complement or reverse-complement of DNA sequences.
#' @param seqs An object of class \code{DNAbin}.
#' @param i Logical or numeric index to indicate a subset of sequences to
#'   manipulate. If not given the entire set of sequences will be manipulated.
#' @param complement Logical, indicating if sequences will be turned into their
#'   complement.
#' @param reverse Logical, indication if sequences will be reversed.
#' @param delete.gaps Logical, indicating if gap character will be removed prior
#'   to sequence manipulation.
#' @return An object of the same class and dimension as \code{seqs}.
#' @examples
#' ## A minimal sequence alignment:
#' x <- list(
#'   seqA = c("a", "a", "c", "c", "g", "t"),
#'   seqB = c("n", "-", "r", "y", "g", "t"))
#' x <- as.DNAbin(x)
#' 
#' ## Three choices of manipulation:
#' as.character(x)
#' as.character(rc(x))                      ## reverse-complement
#' as.character(rc(x, complement = FALSE))  ## only reverse
#' as.character(rc(x, reverse = FALSE))     ## only complement
#' 
#' ## You can remove gaps:
#' as.character(rc(x, delete.gaps = TRUE))  ## gaps/indels removed
#'
#' @export

rc <- function(seqs, i, complement = TRUE, reverse = TRUE, delete.gaps = FALSE){
  
  if (delete.gaps) seqs <- del.gaps(seqs)
  if (missing(i)) i <- 1:length(seqs)
  
  ## IUPAC ambiguity code
  ## --------------------
  iupac <- c(n = 240, "?" = 2, "-" = 4,
             a = 136, c = 40, g = 72, t = 24, 
             r = 192, y = 48, s = 96, w = 144, k = 80, m = 160, 
             b = 112, d = 208, h = 176, v = 224)
  iupac_names <- names(iupac)
  iupac <- as.raw(iupac)
  names(iupac) <- iupac_names
  
  core <- function(seq, complement, reverse){
    
    ## Get complement of sequence
    ## --------------------------
    if (complement){
      seq_rc <- vector(mode = "raw", length = length(seq))
      seq_rc[seq == iupac[["a"]]] <- iupac[["t"]]
      seq_rc[seq == iupac[["t"]]] <- iupac[["a"]]
      seq_rc[seq == iupac[["c"]]] <- iupac[["g"]]
      seq_rc[seq == iupac[["g"]]] <- iupac[["c"]]
      seq_rc[seq == iupac[["r"]]] <- iupac[["y"]] # R = A|G
      seq_rc[seq == iupac[["y"]]] <- iupac[["r"]] # Y = C|T
      seq_rc[seq == iupac[["w"]]] <- iupac[["s"]] # W = A|T
      seq_rc[seq == iupac[["s"]]] <- iupac[["w"]] # S = G|C
      seq_rc[seq == iupac[["m"]]] <- iupac[["k"]] # M = A|C
      seq_rc[seq == iupac[["k"]]] <- iupac[["m"]] # K = G|T
      seq_rc[seq == iupac[["h"]]] <- iupac[["b"]] # H = A|T|C != G
      seq_rc[seq == iupac[["b"]]] <- iupac[["h"]] # B = G|C|T != A
      seq_rc[seq == iupac[["v"]]] <- iupac[["d"]] # V = G|A|C != T
      seq_rc[seq == iupac[["d"]]] <- iupac[["v"]] # D = A|G|T != C
      seq_rc[seq == iupac[["n"]]] <- iupac[["n"]]
      seq_rc[seq == iupac[["-"]]] <- iupac[["-"]]
      if (any(seq_rc %in% as.raw(00))) {
        stop("implement ambiguity codes")
      } 
    } else {
      seq_rc <- as.raw(seq)
    }

    ## Get reverse of sequence
    ## --------------------------
    if (reverse){
      seq_rc <- rev(seq_rc)
    } 
    
    return(seq_rc)
  }
  seqs[i] <- lapply(seqs[i], core, complement = complement, reverse = reverse)
  class(seqs) <- "DNAbin"
  seqs
}
