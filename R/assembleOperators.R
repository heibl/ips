## This code is part of the ips package
## © C. Heibl 2014 (last update 2019-12-12)

assembleOperators <- function(id, link.clocks, link.trees){
  
  o <- list()
  for (i in seq_along(id)){
    
    if (link.trees){
      
      ## *.t operators for unlinked trees
      ## --------------------------------
      if (i == 1){
        ## *.t operators for unlinked trees
        ## --------------------------------
        o <- c(o, list(xmlNode("operator", 
                               attrs = c(id = "YuleBirthRateScaler.t:trees",
                                         parameter = "@birthRate.t:trees",
                                         scaleFactor = "0.75",
                                         spec = "ScaleOperator",
                                         weight = "3.0"))))
        
        
        
        o <- c(o, list(xmlNode("operator", 
                            attrs = c(id = "treeScaler.t:trees",
                                      scaleFactor = "0.5",
                                      spec = "ScaleOperator",
                                      tree = "@Tree.t:trees",
                                      weight = "3.0"))))
        
        
        o <- c(o, list(xmlNode("operator", 
                            attrs = c(id = "treeRootScaler.t:trees",
                                      rootOnly = "true",
                                      scaleFactor = "0.5",
                                      spec = "ScaleOperator",
                                      tree = "@Tree.t:trees",
                                      weight = "3.0"))))
        
        o <- c(o, list(xmlNode("operator", 
                            attrs = c(id = "UniformOperator.t:trees", 
                                      spec = "Uniform",
                                      tree = "@Tree.t:trees",
                                      weight = "30.0"))))
        
        o <- c(o, list(xmlNode("operator", 
                            attrs = c(id = "SubtreeSlide.t:trees",
                                      spec = "SubtreeSlide",
                                      tree = "@Tree.t:trees",
                                      weight = "15.0"))))
        
        
        o <- c(o, list(xmlNode("operator", 
                            attrs = c(id = "narrow.t:trees",
                                      spec = "Exchange",
                                      tree = "@Tree.t:trees",
                                      weight = "15.0"))))
        
        
        o <- c(o, list(xmlNode("operator", 
                            attrs = c(id = "wide.t:trees",
                                      isNarrow = "false",
                                      spec = "Exchange",
                                      tree = "@Tree.t:trees",
                                      weight = "3.0"))))
        
        
        o <- c(o, list(xmlNode("operator", 
                            attrs = c(id = "WilsonBalding.t:trees",
                                      spec = "WilsonBalding",
                                      tree = "@Tree.t:trees",
                                      weight = "3.0"))))
      } else {
        ## Do nothing
      }
      
    } else {
      
      ## *.t operators for unlinked trees
      ## --------------------------------
      
      o <- c(o, list(xmlNode("operator", 
                          attrs = c(id = paste0("YuleBirthRateScaler.t:", id[i]),
                                    parameter = paste0("@birthRate.t:", id[i]),
                                    scaleFactor = "0.75",
                                    spec = "ScaleOperator",
                                    weight = "3.0"))))
      
      
      
      o <- c(o, list(xmlNode("operator", 
                          attrs = c(id = paste0("treeScaler.t:", id[i]),
                                    scaleFactor = "0.5",
                                    spec = "ScaleOperator",
                                    tree = paste0("@Tree.t:", id[i]),
                                    weight = "3.0"))))
      
      
      o <- c(o, list(xmlNode("operator", 
                          attrs = c(id = paste0("treeRootScaler.t:", id[i]),
                                    rootOnly = "true",
                                    scaleFactor = "0.5",
                                    spec = "ScaleOperator",
                                    tree = paste0("@Tree.t:", id[i]),
                                    weight = "3.0"))))
      
      o <- c(o, list(xmlNode("operator", 
                          attrs = c(id = paste0("UniformOperator.t:", id[i]), 
                                    spec = "Uniform",
                                    tree = paste0("@Tree.t:", id[i]),
                                    weight = "30.0"))))
      
      o <- c(o, list(xmlNode("operator", 
                          attrs = c(id = paste0("SubtreeSlide.t:", id[i]),
                                    spec = "SubtreeSlide",
                                    tree = paste0("@Tree.t:", id[i]),
                                    weight = "15.0"))))
      
      
      o <- c(o, list(xmlNode("operator", 
                          attrs = c(id = paste0("narrow.t:", id[i]),
                                    spec = "Exchange",
                                    tree = paste0("@Tree.t:", id[i]),
                                    weight = "15.0"))))
      
      
      o <- c(o, list(xmlNode("operator", 
                          attrs = c(id = paste0("wide.t:", id[i]),
                                    isNarrow = "false",
                                    spec = "Exchange",
                                    tree = paste0("@Tree.t:", id[i]),
                                    weight = "3.0"))))
      
      
      o <- c(o, list(xmlNode("operator", 
                          attrs = c(id = paste0("WilsonBalding.t:", id[i]),
                                    spec = "WilsonBalding",
                                    tree = paste0("@Tree.t:", id[i]),
                                    weight = "3.0"))))
    }
    
    if (!link.clocks & i > 1){
  
      o <- c(o, list(xmlNode("operator", 
                          attrs = c(id = paste0("StrictClockRateScaler.c:", id[i]),
                                    parameter = paste0("@clockRate.c:", id[i]),
                                    scaleFactor = "0.75",
                                    spec = "ScaleOperator",
                                    weight = "3.0"))))
      
      Tree.t <- ifelse(link.trees, "Tree.t:trees", paste0("Tree.t:", id[i]))
      o <- c(o, list(xmlNode("operator", 
                          attrs = c(id = paste0("strictClockUpDownOperator.c:", id[i]),
                                    scaleFactor = "0.75",
                                    spec = "UpDownOperator",
                                    weight = "3.0"),
                          .children = list(xmlNode("parameter",
                                                   attrs = c(idref = paste0("clockRate.c:", id[i]),
                                                             name = "up")),
                                           xmlNode("tree",
                                                   attrs = c(idref = Tree.t,
                                                             name = "down"))))))
    }
  }
  o
}

