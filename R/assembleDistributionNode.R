## This code is part of the ips package
## © C. Heibl 2014 (last update 2025-08-30)

#' @importFrom XML addChildren xmlAttrs<-
#' @export

assembleDistributionNode <- function(x){
  
  xmlNode("distribution", 
          attrs = c(id = "posterior",
                    spec = "util.CompoundDistribution"),
          .children = list(assemblePriorNode(x), 
                           assembleLikelihoodNode(x)))
}