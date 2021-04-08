## This code is part of the ips package
## © C. Heibl 2014 (last update 2020-02-14)

#' @export

assembleInitNode <- function(x){
  
  init <- list()
  for (i in seq_along(x$partition)){
    if (x$link.trees){
      if (i == 1){
        randomPopSize.t <- xmlNode("parameter", "1.0", 
                                   attrs = c(id = paste0("randomPopSize.t:", x$partition[i]),
                                             spec = "parameter.RealParameter",
                                             name = "popSize"))
        ConstantPopulation0.t <- xmlNode("populationModel",
                                         attrs = c(id = paste0("ConstantPopulation0.t:", x$partition[i]),
                                                   spec = "ConstantPopulation"),
                                         .children = list(randomPopSize.t))
        init <- c(init, list(
          xmlNode("init", 
                  attrs = c(id = paste0("RandomTree.t:", x$partition[i]),
                            spec = "beast.evolution.tree.RandomTree",
                            estimate = "false",
                            initial = paste0("@Tree.t:", x$partition[i]),
                            taxa = paste0("@", x$partition[i])),
                  .children = list(ConstantPopulation0.t))))
      } else {
        ## do nothing
      }
    } else {
      randomPopSize.t <- xmlNode("parameter", "1.0", 
                                 attrs = c(id = paste0("randomPopSize.t:", x$partition[i]),
                                           spec = "parameter.RealParameter",
                                           name = "popSize"))
      ConstantPopulation0.t <- xmlNode("populationModel",
                                       attrs = c(id = paste0("ConstantPopulation0.t:", x$partition[i]),
                                                 spec = "ConstantPopulation"),
                                       .children = list(randomPopSize.t))
      init <- c(init, list(
        xmlNode("init", 
                attrs = c(id = paste0("RandomTree.t:", x$partition[i]),
                          spec = "beast.evolution.tree.RandomTree",
                          estimate = "false",
                          initial = paste0("@Tree.t:", x$partition[i]),
                          taxa = paste0("@", x$partition[i])),
                .children = list(ConstantPopulation0.t))))
    }
  }
  init
}
