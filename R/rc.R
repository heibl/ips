## This code is part of the ips package
## © C. Heibl 2018 (last update 2019-10-15)

#' @export

rc <- function(seqs, i){
  
  seqs <- del.gaps(seqs)
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
  
  core <- function(seq){
    rc_seq <- vector(mode = "raw", length = length(seq))
    rc_seq[seq == iupac[["a"]]] <- iupac[["g"]]
    rc_seq[seq == iupac[["g"]]] <- iupac[["a"]]
    rc_seq[seq == iupac[["c"]]] <- iupac[["t"]]
    rc_seq[seq == iupac[["t"]]] <- iupac[["c"]]
    rc_seq[seq == iupac[["r"]]] <- iupac[["y"]] # R = A|G
    rc_seq[seq == iupac[["y"]]] <- iupac[["r"]] # Y = C|T
    rc_seq[seq == iupac[["w"]]] <- iupac[["s"]] # W = A|T
    rc_seq[seq == iupac[["s"]]] <- iupac[["w"]] # S = G|C
    rc_seq[seq == iupac[["m"]]] <- iupac[["k"]] # M = A|C
    rc_seq[seq == iupac[["k"]]] <- iupac[["m"]] # K = G|T
    rc_seq[seq == iupac[["h"]]] <- iupac[["b"]] # H = A|T|C != G
    rc_seq[seq == iupac[["b"]]] <- iupac[["h"]] # B = G|C|T != A
    rc_seq[seq == iupac[["v"]]] <- iupac[["d"]] # V = G|A|C != T
    rc_seq[seq == iupac[["d"]]] <- iupac[["v"]] # D = A|G|T != C
    rc_seq[seq == iupac[["n"]]] <- iupac[["n"]]
    if (any(rc_seq %in% as.raw(00))) {
      stop("implement ambiguity codes")
    }
    rev(rc_seq)
  }
  seqs[i] <- lapply(seqs[i], core)
  class(seqs) <- "DNAbin"
  seqs
}
