## This code is part of the ips package
## © C. Heibl 2014 (last update 2015-04-05)

assembleDistributionNode <- function(id){
  
  prior <- xmlNode("distribution", 
                   attrs = c(id = "prior",
                             spec = "util.CompoundDistribution"))
  u <- xmlNode("Uniform",
               attrs = c(id = "Uniform.0",
                         name ="distr",
                         upper = "Infinity"))
  u.step <- 1
  
  ## node likelihood
  likelihood <- xmlNode(
    "distribution", 
    attrs = c(id = "likelihood",
              spec = "util.CompoundDistribution"))
  
  for ( i in seq_along(id) ){
    
    ## add children to prior
    ## ---------------------
    
    d <- xmlNode("distribution", 
                 attrs = c(birthDiffRate = paste("@birthRate.t:", id[i], sep = ""),
                           id = paste("YuleModel.t:", id[i], sep = ""),
                           spec = "beast.evolution.speciation.YuleModel",
                           tree = paste("@Tree.t:", id[i], sep = "")))
    prior <- addChildren(prior, kids = list(d))
    
    if ( i > 1 ) {
      xmlAttrs(u) <- c(id = paste("Uniform.0", u.step, sep = ""))
      u.step <- u.step + 1
      
      p <- xmlNode("prior",
                   attrs = c(id = paste("ClockPrior.c:", id[i], sep = ""),
                             name = "distribution",
                             x = paste("@clockRate.c:", id[i], sep = "")),
                   .children = list(u))
      prior <- addChildren(prior, kids = list(p))
      
      xmlAttrs(u) <- c(id = paste("Uniform.0", u.step, sep = ""))
      u.step <- u.step + 1
    }
    p <- xmlNode("prior",
                 attrs = c(id = paste("YuleBirthRatePrior.t:", id[i], sep = ""),
                           name = "distribution",
                           x = paste("@birthRate.t:", id[i], sep = "")),
                 .children = list(u))
    prior <- addChildren(prior, kids = list(p))
    
    ## add children to node likelihood
    ## -------------------------------
    siteModel <- xmlNode(
      "siteModel", 
      attrs = c(id = paste("SiteModel.s:", id[i], sep = ""),
                spec = "SiteModel"),
      .children = list(xmlNode("parameter",  "1.0",
                               attrs = c(estimate = "false",
                                         id = paste("mutationRate.s:", id[i], sep = ""),
                                         name = "mutationRate")),
                       xmlNode("parameter", "1.0",
                               attrs = c(estimate = "false",
                                         id = paste("gammaShape.s:", id[i], sep = ""),
                                         name = "shape")),
                       xmlNode("parameter", "0.0", 
                               attrs = c(estimate = "false",
                                         id = paste("proportionInvariant.s:", id[i], sep = ""),
                                         lower = "0.0",
                                         name = "proportionInvariant",
                                         upper = "1.0")),
                       xmlNode("substModel", attrs = c(id = paste("JC69.s:", id[i], sep = ""),
                                                       spec="JukesCantor"))))
    if ( i == 1){
      branchRateModel <- xmlNode(
        "branchRateModel",
        attrs = c(id = paste("StrictClock.c:", id[i], sep = ""),
                  spec = "beast.evolution.branchratemodel.StrictClockModel"),
        xmlNode("parameter", "1.0", 
                attrs = c(estimate = "false",
                          id = paste("clockRate.c:", id[i], sep = ""),
                          name = "clock.rate")))
    } else {
      branchRateModel <- xmlNode(
        "branchRateModel",
        attrs = c(clock.rate = paste("@clockRate.c:", id[i], sep = ""),
                  id = paste("StrictClock.c:", id[i], sep = ""),
                  spec = "beast.evolution.branchratemodel.StrictClockModel"))
    }
    
    distribution <- xmlNode(
      "distribution",
      attrs = c(data = paste("@", id[i], sep = ""),
                id = paste("treeLikelihood.", id[i], sep = ""),
                spec = "TreeLikelihood",
                tree = paste("@Tree.t:", id[i], sep = "")),
      .children = list(siteModel, branchRateModel))
    likelihood <- addChildren(likelihood, kids = list(distribution))
  }
  
  
  
  #   if ( !missing(taxonset) ) {
  #     ts <- lapply(taxonset, rbeauti.taxonset, id = id)
  #     prior.distribution <- addChildren(prior.distribution, kids = ts)
  #   }
  
  
  xmlNode("distribution", 
          attrs = c(id = "posterior",
                    spec = "util.CompoundDistribution"),
          .children = list(prior, likelihood))
  
}