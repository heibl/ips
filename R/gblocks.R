## This code is part of the ips package
## © C. Heibl 2014 (last update 2016-11-08)

#' @export

gblocks <- function(x, b1 = .5, b2 = b1, b3 = ncol(x), 
                    b4 = 2, b5 = "a", exec){
  
  if ( inherits(x, "alignment") ) stop("cannot handle class 'alignment'")
  if ( inherits(x, "list") ) stop("cannot handle unaligned sequences")
  
  ## check parameters:
  ## -----------------
  if ( b1 < .5 | b1 > 1 ) stop ("b1 not in [0.5, 1]")
  if ( b2 < b1 | b2 > 1 ) stop ("b2 not in [b1, 1]")
  if ( b3 < 0 | b4 > ncol(x) ) stop ("b3 not in [0, ", ncol(x), "]")
  if ( b4 < 2 | b4 > ncol(x) ) stop ("b4 not in [2, ", ncol(x), "]")
  b5 <- match.arg(b5, c("a", "h", "n"))
  
  if ( !file.exists(exec) ) 
    stop("executable '", exec, "' does not exist", sep = "")
  
  
  ntax <- nrow(x)
  b1 <- floor(ntax * b1) + 1
  b2 <- floor(ntax * b2) + 1
  
  cat("\n--- executing Gblocks ---")
  cat("\nminimum number of sequences for a conserved position :", b1)
  cat("\nminimum number of sequences for a flank position     :", b2)
  cat("\nmaximum number of contiguous nonconserved positions  :", b3)
  cat("\nminimum length of a block                            :", b4)
  cat("\nallowed gap positions                                :", b5)
  
  write.fas(x, "R2GBLOCK.fas")
  system(paste(exec, " R2GBLOCK.fas -t=d", 
               " -b1=", b1,
               " -b2=", b2,
               " -b3=", b3,
               " -b4=", b4,
               " -b5=", b5, 
               sep = "")
         )
  out <- read.fas("R2GBLOCK.fas-gb")
  unlink(list.files(pattern = "R2GBLOCK"))
  out
}