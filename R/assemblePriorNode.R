## This code is part of the ips package
## © C. Heibl 2020 (last update 2020-03-13)

assemblePriorNode <- function(x){
  
  # ## Loop over partitions to add childen to <prior>
  # ## ----------------------------------------------
  # for (i in seq_along(id)){
  #   prior <- assemblePriorNode(id, i, link.trees, prior)
  # }
  i <- 1
  
  # ## Define taxonsets
  # ## ----------------
  # if (!missing(taxonset)) {
  #   
  #   ## Loop over taxonsets
  #   ## -------------------
  #   for (j in seq_along(taxonset)){
  #     
  #     taxonSetPrior <- assembleTaxonSetPrior(taxonset = taxonset[[j]],
  #                                            setname = names(taxonset)[j],
  #                                            id = id, i = i, link.trees = link.trees)
  #     if (!is.null(taxonSetPrior)){
  #       prior <- addChildren(prior, kids = list(taxonSetPrior))
  #     }
  #   }
  # }
  
  prior <- xmlNode("distribution", 
                   attrs = c(id = "prior",
                             spec = "util.CompoundDistribution"))
  
  if (x$tree == "Yule"){
    kids <- list(
      xmlNode("distribution", 
              attrs = c(id = paste0("YuleModel.t:", x$partition[i]),
                        spec = "beast.evolution.speciation.YuleModel",
                        birthDiffRate = paste0("@birthRate.t:", x$partition[i]),
                        tree = paste0("@Tree.t:", x$partition[i]))),
      xmlNode("prior",
              attrs = c(id = paste0("YuleBirthRatePrior.t:", x$partition[i]),
                        name = "distribution",
                        x = paste0("@birthRate.t:", x$partition[i])),
              .children = uniformPrior(x))
    )
  }
  
  if (x$tree == "Calibrated Yule"){
    kids <- list(
      xmlNode("distribution", 
              attrs = c(birthRate = paste0("@birthRateY.t:", x$partition[i]),
                        id = paste0("CalibratedYuleModel.t:", x$partition[i]),
                        spec = "beast.evolution.speciation.CalibratedYuleModel",
                        tree = paste0("@Tree.t:", x$partition[i]))),
      xmlNode("prior",
              attrs = c(id = paste0("CalibratedYuleBirthRatePrior.t:", x$partition[i]),
                        name = "distribution",
                        x = paste0("@birthRateY.t:", x$partition[i])),
              .children = uniformPrior(x, upper = "1000.0"))
    )
  }
  
  if (x$tree == "Birth Death"){
    kids <- list(
      xmlNode("distribution", 
              attrs = c(birthDiffRate = paste0("@BDBirthRate.t:", x$partition[i]),
                        id = paste0("BirthDeath.t:", x$partition[i]),
                        relativeDeathRate = paste0("@BDDeathRate.t:", x$partition[i]),
                        spec = "beast.evolution.speciation.BirthDeathGernhard08Model",
                        tree = paste0("@Tree.t:", x$partition[i]))),
      xmlNode("prior",
              attrs = c(id = paste0("BirthRatePrior.t:", x$partition[i]),
                        name = "distribution",
                        x = paste0("@BDBirthRate.t:", x$partition[i])),
              .children = uniformPrior(x, upper = "1000.0")),
      xmlNode("prior",
              attrs = c(id = paste0("DeathRatePrior.t:", x$partition[i]),
                        name = "distribution",
                        x = paste0("@BDDeathRate.t:", x$partition[i])),
              .children = uniformPrior(x, upper = NULL))
    )
  }
  
  if (x$tree == "Fossilized Birth Death"){
    kids <- list(
      xmlNode("distribution", 
              attrs = c(id = paste0("FBD.t:", x$partition[i]),
                        spec = "beast.evolution.speciation.SABirthDeathModel",
                        conditionOnRhoSampling = "true",
                        diversificationRate = paste0("@diversificationRateFBD.t:", x$partition[i]),
                        origin = paste0("@originFBD.t:", x$partition[i]),
                        samplingProportion = paste0("@samplingProportionFBD.t:", x$partition[i]),
                        tree = paste0("@Tree.t:", x$partition[i]),
                        turnover = paste0("@turnoverFBD.t:", x$partition[i])),
              .children = lapply(c("rFBD.t", "rhoFBD.t"), parameter, i = i, x = x)),
      xmlNode("prior",
              attrs = c(id = paste0("diversificationRatePriorFBD.t:", x$partition[i]),
                        name = "distribution",
                        x = paste0("@diversificationRateFBD.t:", x$partition[i])),
              .children = uniformPrior(x)),
      xmlNode("prior",
              attrs = c(id = paste0("originPriorFBD.t:", x$partition[i]),
                        name = "distribution",
                        x = paste0("@originFBD.t:", x$partition[i])),
              .children = uniformPrior(x)),
      xmlNode("prior",
              attrs = c(id = paste0("samplingProportionPriorFBD.t:", x$partition[i]),
                        name = "distribution",
                        x = paste0("@samplingProportionFBD.t:", x$partition[i])),
              .children = uniformPrior(x, upper = NULL)),
      xmlNode("prior",
              attrs = c(id = paste0("turnoverPriorFBD.t:", x$partition[i]),
                        name = "distribution",
                        x = paste0("@turnoverFBD.t:", x$partition[i])),
              .children = uniformPrior(x, upper = NULL)))
    
    if (x$clock == "Strict Clock"){
      
      kids <- c(kids, list(xmlNode("prior",
                                   attrs = c(id = paste0("ClockPrior.c:", x$partition[i]),
                                             name = "distribution",
                                             x = paste0("@clockRate.c:", x$partition[i])),
                                   .children = uniformPrior(x))))
    }
    if (x$clock == "Relaxed Clock Exponential"){
      
      UCMeanRatePrior.c <- xmlNode("prior",
                                   attrs = c(id = paste0("UCMeanRatePrior.c:", x$partition[i]),
                                             name = "distribution",
                                             x = paste0("@ucedMean.c:", x$partition[i])),
                                   .children = uniformPrior(x))
      kids <- c(kids, list(UCMeanRatePrior.c))
      
    }
    
  }
  
  ## RELAXED CLOCK LOG NORMAL  
  ## ------------------------
  if (x$clock == "Relaxed Clock Log Normal"){
    
    ## Default: 'ucldStdevPrior.c'
    clock.kids <- list(xmlNode("prior",
                               attrs = c(id = paste0("ucldStdevPrior.c:", x$partition[i]),
                                         name = "distribution",
                                         x = paste0("@ucldStdev.c:", x$partition[i])),
                               .children = list(gammaPrior(x, "0.5396", "0.3819"))))
    
    ## Have to add 'MeanRatePrior.c', when tip dating                                                                                                             .children = list(0.3819)))))))
    if (!is.na(x$tip.dates)){
      clock.kids <- c(list(xmlNode("prior", attrs = c(id = paste0("MeanRatePrior.c:", x$partition[i]),
                                                      name = "distribution",
                                                      x = paste0("@ucldMean.c:", x$partition[i])),
                                   .children = uniformPrior(x))),
                      clock.kids)
    }
    kids <- c(kids, clock.kids)
  }
  
  ## RANDOM LOCAL CLOCK  
  ## ------------------------
  if (x$clock == "Random Local Clock"){
    clock.kids <- list(xmlNode("prior",
                               attrs = c(id = paste0("RRateChangesPrior.c:", x$partition[i]),
                                         name = "distribution"),
                               .children = list(
                                 xmlNode("x", attrs = c(id = paste0("RRateChanges.c:", x$partition[i]),
                                                      spec = "util.Sum"),
                                         value = xmlNode("arg", attrs = c(idref = paste0("Indicators.c:", x$partition[i])))),
                                 xmlNode("distr", attrs = c(id = "Poisson.0",
                                                            spec = "beast.math.distributions.Poisson"),
                                         value = xmlNode("parameter", 
                                                         attrs = c(id = realParameter(x),
                                                                   spec = "parameter.RealParameter",
                                                                   estimate = "false",
                                                                   name = "lambda"),
                                                         value = "0.6931471805599453"))
                               )),
                       xmlNode("prior",
                               attrs = c(id = paste0("RRatesPrior.c:", x$partition[i]),
                                         name = "distribution",
                                         x = paste0("@clockrates.c:", x$partition[i])),
                               .children = list(gammaPrior(x, "0.5396", "0.3819"))))
    kids <- c(kids, clock.kids)
  }
  
  if (x$subst.model == "HKY"){
    
    FrequenciesPrior.s <- xmlNode("prior", 
                                  attrs = c(id = paste0("FrequenciesPrior.s:", x$partition[i]),
                                            name = "distribution",
                                            x = paste0("@freqParameter.s:", x$partition[i])),
                                  .children = uniformPrior(x, upper = NULL))
    
    KappaPrior.s <- list(
      xmlNode("parameter", 
              attrs = c(estimate = "false",
                        id = "RealParameter.17",
                        name = "M",
                        spec = "parameter.RealParameter"),
              value = "1.0"),
      xmlNode("parameter", 
              attrs = c(estimate = "false",
                        id = "RealParameter.18",
                        name = "S",
                        spec = "parameter.RealParameter"),
              value = "1.25"))
    KappaPrior.s <- xmlNode("LogNormal",
                            attrs = c(id = "LogNormalDistributionModel.2", name = "distr"),
                            .children = KappaPrior.s)
    KappaPrior.s <- xmlNode("prior",
                            attrs = c(id = paste0("KappaPrior.s:", x$partition[i]),
                                      name = "distribution",
                                      x = paste0("@kappa.s:", x$partition[i])),
                            value = KappaPrior.s)
    
    kids <- c(kids, list(FrequenciesPrior.s, KappaPrior.s))
  }
  
  if (x$subst.model == "TN93"){
    
    FrequenciesPrior.s <- xmlNode("prior", 
                                  attrs = c(id = paste0("FrequenciesPrior.s:", x$partition[i]),
                                            name = "distribution",
                                            x = paste0("@freqParameter.s:", x$partition[i])),
                                  .children = uniformPrior(x, upper = NULL))
    
    kappa1Prior.s <- xmlNode("prior",
                            attrs = c(id = paste0("kappa1Prior.s:", x$partition[i]),
                                      name = "distribution",
                                      x = paste0("@kappa1.s:", x$partition[i])),
                            value = LogNormalPrior(x, 1.0, 1.25))
    
    kappa2Prior.s <- xmlNode("prior",
                             attrs = c(id = paste0("kappa2Prior.s:", x$partition[i]),
                                       name = "distribution",
                                       x = paste0("@kappa2.s:", x$partition[i])),
                             value = LogNormalPrior(x, 1.0, 1.25))
    
    
    kids <- c(kids, list(FrequenciesPrior.s, kappa1Prior.s, kappa2Prior.s))
  }
  
  if (x$subst.model == "GTR"){
    
    FrequenciesPrior.s <- xmlNode("prior", 
                                  attrs = c(id = paste0("FrequenciesPrior.s:", x$partition[i]),
                                            name = "distribution",
                                            x = paste0("@freqParameter.s:", x$partition[i])),
                                  .children = uniformPrior(x, upper = NULL))
    
    RateACPrior.s <- xmlNode("prior", 
                             attrs = c(id = paste0("RateACPrior.s:", x$partition[i]),
                                       name = "distribution",
                                       x = paste0("@rateAC.s:", x$partition[i])),
                             value = gammaPrior(x, "0.05", "10.0"))
    RateAGPrior.s <- xmlNode("prior", 
                             attrs = c(id = paste0("RateAGPrior.s:", x$partition[i]),
                                       name = "distribution",
                                       x = paste0("@rateAG.s:", x$partition[i])),
                             value = gammaPrior(x, "0.05", "20.0")) ## RLC beta of AC and AG exchanged in BEAUTI
    RateATPrior.s <- xmlNode("prior", 
                             attrs = c(id = paste0("RateATPrior.s:", x$partition[i]),
                                       name = "distribution",
                                       x = paste0("@rateAT.s:", x$partition[i])),
                             value = gammaPrior(x, "0.05", "10.0"))
    RateCGPrior.s <- xmlNode("prior", 
                             attrs = c(id = paste0("RateCGPrior.s:", x$partition[i]),
                                       name = "distribution",
                                       x = paste0("@rateCG.s:", x$partition[i])),
                             value = gammaPrior(x, "0.05", "10.0"))
    RateGTPrior.s <- xmlNode("prior", 
                             attrs = c(id = paste0("RateGTPrior.s:", x$partition[i]),
                                       name = "distribution",
                                       x = paste0("@rateGT.s:", x$partition[i])),
                             value = gammaPrior(x, "0.05", "10.0"))
    
    kids <- c(kids, list(FrequenciesPrior.s, RateACPrior.s, RateAGPrior.s, 
                         RateATPrior.s, RateCGPrior.s, RateGTPrior.s))
  }
  
  prior <- addChildren(prior, kids = kids)
  prior
}


