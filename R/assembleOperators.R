## This code is part of the ips package
## © C. Heibl 2014 (last update 2015-04-05)

assembleOperators <- function(id){
  
  ops <- vector(mode = "list")
  for ( i in seq_along(id) ){
    this.ops <- list(
             xmlNode("operator", 
                     attrs = c(id = paste("YuleBirthRateScaler.t:", id[i], sep = ""),
                               parameter = paste("@birthRate.t:", id[i], sep = ""),
                               scaleFactor = "0.75",
                               spec = "ScaleOperator",
                               weight = "3.0")),
             xmlNode("operator", 
                     attrs = c(id = paste("treeScaler.t:", id[i], sep = ""),
                               scaleFactor = "0.5",
                               spec = "ScaleOperator",
                               tree = paste("@Tree.t:", id[i], sep = ""),
                               weight = "3.0")),
             xmlNode("operator", 
                     attrs = c(id = paste("treeRootScaler.t:", id[i], sep = ""),
                               rootOnly = "true",
                               scaleFactor = "0.5",
                               spec = "ScaleOperator",
                               tree = paste("@Tree.t:", id[i], sep = ""),
                               weight = "3.0")),
             xmlNode("operator", 
                     attrs = c(id = paste("UniformOperator.t:", id[i], sep = ""), 
                               spec = "Uniform",
                               tree = paste("@Tree.t:", id[i], sep = ""),
                               weight = "30.0")),
             xmlNode("operator", 
                     attrs = c(id = paste("SubtreeSlide.t:", id[i], sep = ""),
                               spec = "SubtreeSlide",
                               tree = paste("@Tree.t:", id[i], sep = ""),
                               weight = "15.0")),
             xmlNode("operator", 
                     attrs = c(id = paste("narrow.t:", id[i], sep = ""),
                               spec = "Exchange",
                               tree = paste("@Tree.t:", id[i], sep = ""),
                               weight = "15.0")),
             xmlNode("operator", 
                     attrs = c(id = paste("wide.t:", id[i], sep = ""),
                               isNarrow = "false",
                               spec = "Exchange",
                               tree = paste("@Tree.t:", id[i], sep = ""),
                               weight = "3.0")),
             xmlNode("operator", 
                     attrs = c(id = paste("WilsonBalding.t:", id[i], sep = ""),
                               spec = "WilsonBalding",
                               tree = paste("@Tree.t:", id[i], sep = ""),
                               weight = "3.0"))
             
    )
    if ( i == 1 ){
      ops <- c(ops, this.ops)
    } else {
      ops <- c(ops, 
               list(xmlNode("operator", 
                            attrs = c(id = paste("StrictClockRateScaler.c:", id[i], sep = ""),
                                      parameter = paste("@clockRate.c:", id[i], sep = ""),
                                      scaleFactor = "0.75",
                                      spec = "ScaleOperator",
                                      weight = "3.0"))),
               this.ops,
               list(xmlNode("operator", 
                            attrs = c(id = paste("strictClockUpDownOperator.c:", id[i], sep = ""),
                                      scaleFactor = "0.75",
                                      spec = "UpDownOperator",
                                      weight = "3.0"),
                            .children = list(xmlNode("parameter",
                                                     attrs = c(idref = paste("clockRate.c:", id[i], sep = ""),
                                                               name = "up")),
                                             xmlNode("tree",
                                                     attrs = c(idref = paste("Tree.t:", id[i], sep = ""),
                                                               name = "down")))))
               )
    }
    
  }
  ops
}
  
  