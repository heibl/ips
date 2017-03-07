## This code is part of the ips package
## © C. Heibl 2014 (last update 2015-04-04)

assembleInitNode <- function(id){
  
  parameter <- xmlNode("parameter", "1.0", 
                       attrs = c(id = paste("randomPopSize.t:", id, sep = ""),
                                 name = "popSize"))
  populationModel <- xmlNode("populationModel",
                             attrs = c(id = paste("ConstantPopulation0.t:", id, sep = ""),
                                       spec = "ConstantPopulation"),
                             .children = list(parameter))
  xmlNode("init", 
          attrs = c(estimate = "false",
                    id = paste("RandomTree.t:", id, sep = ""),
                    initial = paste("@Tree.t:", id, sep = ""),
                    spec = "beast.evolution.tree.RandomTree",
                    taxa = paste("@", id, sep = "")),
          .children = list(populationModel))
  
}
