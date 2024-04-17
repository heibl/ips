## This code is part of the ips package
## © C. Heibl 2014 (last update 2015-04-05)

setClock  <- function(xml, clock){
  
  clock <- match.arg(clock, c("Strict Clock",
                              "Relaxed Clock Exponential",
                              "Relaxed Clock Log Normal",
                              "Random Local Clock"))
  
  
  if ( clock == "Relaxed Clock Exponential") {
    xml[["run"]][["state"]] <- addChildren(
      xml[["run"]][["state"]],
      kids = list(xmlNode("stateNode", 1,
                          attrs = c(dimension = "8",
                                    id = paste("expRateCategories.c:", "xxx", sep = ""),
                                    spec = "parameter.IntegerParameter"))))
    ## branchRateModel
    ## ---------------
    p <- xmlNode("parameter", "1.0",
                 attrs = c(id="UCExpLambda.c:xxx1",
                           name="mean"))
    p <- xmlNode("Exponential", p,
                 attrs = c(id = "Exponential.c:xxx1",
                           name="distr"))
    branchRateModel <- xmlNode(
      "branchRateModel",
      attrs = c(id = "ExponentialRelaxedClock.c:xxx1", 
                rateCategories = "@expRateCategories.c:xxx1",
                spec="beast.evolution.branchratemodel.UCRelaxedClockModel",
                tree="@Tree.t:xxx1"),
      .children = list(p,
                  xmlNode("parameter", "1.0",
                          attrs = c(estimate = "false",
                                    id = "ucedMean.c:xxx1",
                                    name = "clock.rate"))))
    xml[["run"]][["distribution"]][[2]][[1]][["branchRateModel"]] <- branchRateModel
  
    ## operators
    ## ---------
    o <- xmlNode("operator",
                 attrs = c(id = "ExpCategoriesRandomWalk.c:xxx1",
                           parameter="@expRateCategories.c:xxx1",
                           spec="IntRandomWalkOperator",
                           weight="10.0",
                           windowSize="1"))
    
    which(names(xmlChildren(xml[["run"]])) == "operator")
  }
  xml
}

