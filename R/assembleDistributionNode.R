## This code is part of the ips package
## © C. Heibl 2014 (last update 2020-02-14)

#' @export

assembleDistributionNode <- function(x){
  
  xmlNode("distribution", 
          attrs = c(id = "posterior",
                    spec = "util.CompoundDistribution"),
          .children = list(assemblePriorNode(x), 
                           assembleLikelihoodNode(x)))
}