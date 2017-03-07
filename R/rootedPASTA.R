## This code is part of the ips package
## © C. Heibl 2016 (last update 2016-11-15)

## WORK UNDER PROGRESS
## - cannot handle polytomies
## - does not incoporate branch length information

#' @title PASTA with a Rooted Guide Tree
#' @description Implements the PASTA (ultra-large multiple sequence alignment) 
#' algorithm of Mirarab, Nguyen, and Warnow (2014) with a rooted guide tree.
#' @param phy An object of class \code{\link{phylo}}.
#' @param k An integer giving the maximum size of the subgroups that will be 
#' aligned with MAFFT.
#' @references 
#' Mirarab, S., N. Nguyen, and T. Warnow. 2014. PASTA: Ultra-large multiple
#' sequence alignment. \emph{RECOMB} \strong{}: 177-191.
#' @export

rootedPASTA <- function(phy, k){
  
}