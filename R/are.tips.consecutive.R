## This code is part of the ips package
## © C. Heibl 2014

are.tips.consecutive <- function(phy){
  canonical <- seq_along(phy$tip.label)
  given <- phy$edge[, 2][phy$edge[, 2] %in% canonical]
  given <- as.integer(given)
  if ( !identical(canonical, given) ) 
    stop("tips are not numbered consecutively.",
         " Type '?fixNodes' for help.")
}