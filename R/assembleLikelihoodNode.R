## This code is part of the ips package
## © C. Heibl 2020 (last update 2020-03-23)

assembleLikelihoodNode <- function(x){
  
  i <- 1
  
  ## Node <likelihood>
  likelihood <- xmlNode("distribution", 
                        attrs = c(id = "likelihood",
                                  spec = "util.CompoundDistribution",
                                  useThreads = "true"))
  
  ## Add children to node <likelihood>
  ## ---------------------------------
  treeLikelihood <- xmlNode("distribution",
                            attrs = c(id = paste0("treeLikelihood.", x$partition[i]),
                                      spec = "ThreadedTreeLikelihood",
                                      data = paste0("@", x$partition[i]),
                                      tree = paste0("@Tree.t:", x$partition[i])))
  
  
  # if (link.clocks & i > 1){
  #   # a[3] <- "@StrictClock.c:clocks"
  #   # names(a)[3] <- "branchRateModel"
  #   a <- c(a[1:3], branchRateModel = "@StrictClock.c:clocks", a[4])
  # } 
  # distribution <- xmlNode("distribution", attrs = a)
  
  ## SiteModel.s
  ## -----------
  if (x$subst.model == "JC69"){
    subst.model <- xmlNode("substModel",
                           attrs = c(id = paste0("JC69.s:", x$partition[i]),
                                     spec = "JukesCantor"))
  }
  if (x$subst.model == "HKY"){
    subst.model <- xmlNode("substModel",
                           attrs = c(id = paste0("hky.s:", x$partition[i]),
                                     kappa = paste0("@kappa.s:", x$partition[i]),
                                     spec = "HKY"),
                           value = xmlNode("frequencies", attrs = c(frequencies = paste0("@freqParameter.s:", x$partition[i]),
                                                                    id = paste0("estimatedFreqs.s:", x$partition[i]),
                                                                    spec = "Frequencies")))
  }
  if (x$subst.model == "TN93"){
    subst.model <- xmlNode("substModel",
                           attrs = c(id = paste0("tn93.s:", x$partition[i]),
                                     kappa1 = paste0("@kappa1.s:", x$partition[i]),
                                     kappa2 = paste0("@kappa2.s:", x$partition[i]),
                                     spec = "TN93"),
                           value = xmlNode("frequencies", attrs = c(frequencies = paste0("@freqParameter.s:", x$partition[i]),
                                                                    id = paste0("estimatedFreqs.s:", x$partition[i]),
                                                                    spec = "Frequencies")))
  }
  if (x$subst.model == "GTR"){
    subst.model <- list(xmlNode("parameter",
                                attrs = c(id = paste0("rateCT.s:", x$partition[i]),
                                          lower = "0.0",
                                          estimate = "false",
                                          name = "rateCT",
                                          spec = "parameter.RealParameter"),
                                value = "1.0"),
                         xmlNode("frequencies", 
                                 attrs = c(frequencies = paste0("@freqParameter.s:", x$partition[i]),
                                                          id = paste0("estimatedFreqs.s:", x$partition[i]),
                                                          spec = "Frequencies")))
    subst.model <- xmlNode("substModel",
                           attrs = c(id = paste0("gtr.s:", x$partition[i]),
                                     rateAC = paste0("@rateAC.s:", x$partition[i]),
                                     rateAG = paste0("@rateAG.s:", x$partition[i]),
                                     rateAT = paste0("@rateAT.s:", x$partition[i]),
                                     rateCG = paste0("@rateCG.s:", x$partition[i]),
                                     rateGT = paste0("@rateGT.s:", x$partition[i]),
                                     spec = "GTR"),
                           .children = subst.model)
  }
  
 
  siteModel <- xmlNode("siteModel", 
                       attrs = c(id = paste0("SiteModel.s:", x$partition[i]),
                                 spec = "SiteModel"))
  params <- c("mutationRate.s", "gammaShape.s", "proportionInvariant.s")
  siteModel <- addChildren(siteModel, 
                           kids = c(lapply(params, parameter, x = x, i = i),
                                    list(subst.model)))
  ## branchRateModel
  ## ---------------
  if (x$clock == "Strict Clock"){
    if (x$tree %in% c("Yule", "Calibrated Yule", "Birth Death")){
      branchRateModel <- xmlNode("branchRateModel",
                                 attrs = c(id = paste0("StrictClock.c:", x$partition[i]),
                                           spec = "beast.evolution.branchratemodel.StrictClockModel"),
                                 .children = list(xmlNode("parameter",
                                                          attrs = c(estimate = "false",
                                                                    id = paste0("clockRate.c:", x$partition[i]),
                                                                    name = "clock.rate",
                                                                    spec = "parameter.RealParameter"),
                                                          value = "1.0"))
      )
    }
    if (x$tree == "Fossilized Birth Death"){
      branchRateModel <- xmlNode("branchRateModel",
                                 attrs = c(clock.rate = paste0("@clockRate.c:", x$partition[i]),
                                           id = paste0("StrictClock.c:", x$partition[i]),
                                           spec = "beast.evolution.branchratemodel.StrictClockModel"))
    }
  } 
  if (x$clock == "Relaxed Clock Exponential") {
    if (x$tree %in% c("Yule", "Calibrated Yule", "Birth Death")){
      branchRateModel <- 
        xmlNode("branchRateModel",
                attrs = c(id = paste0("ExponentialRelaxedClock.c:", x$partition[i]),
                          rateCategories = paste0("@expRateCategories.c:", x$partition[i]),
                          spec = "beast.evolution.branchratemodel.UCRelaxedClockModel",
                          tree = paste0("@Tree.t:", x$partition[i])),
                .children = list(xmlNode("Exponential",
                                         attrs = c(id = paste0("Exponential.c:", x$partition[i]),
                                                   name = "distr"),
                                         value = parameter("UCExpLambda.c", x = x, i = i)),
                                 parameter("ucedMean.c.Yule", x = x, i = i)))
    }
    if (x$tree == "Fossilized Birth Death"){
      branchRateModel <- 
        xmlNode("branchRateModel",
                attrs = c(id = paste0("ExponentialRelaxedClock.c:", x$partition[i]),
                          rateCategories = paste0("@expRateCategories.c:", x$partition[i]),
                          spec = "beast.evolution.branchratemodel.UCRelaxedClockModel",
                          clock.rate = paste0("@ucedMean.c:", x$partition[i]), ## only difference to Yule
                          tree = paste0("@Tree.t:", x$partition[i])),
                .children = list(xmlNode("Exponential",
                                         attrs = c(id = paste0("Exponential.c:", x$partition[i]),
                                                   name = "distr"),
                                         value = parameter("UCExpLambda.c", x = x, i = i))))
    }
  }
  if (x$clock == "Relaxed Clock Log Normal") {
    if (x$tree %in% c("Yule", "Calibrated Yule", "Birth Death")){
      id <- get("counter", envir = x$environment)
      branchRateModel <- 
        xmlNode("branchRateModel",
                attrs = c(id = paste0("RelaxedClock.c:", x$partition[i]),
                          spec = "beast.evolution.branchratemodel.UCRelaxedClockModel",
                          # clock.rate = paste0("@ucldMean.c:", x$partition[i]), ## only tip dates
                          rateCategories = paste0("@rateCategories.c:", x$partition[i]),
                          tree = paste0("@Tree.t:", x$partition[i])),
                .children = list(xmlNode("LogNormal",
                                         attrs = c(id = paste0("LogNormalDistributionModel.c:", x$partition[i]),
                                                   S = paste0("@ucldStdev.c:", x$partition[i]),
                                                   meanInRealSpace = "true",
                                                   name = "distr"),
                                         value = xmlNode("parameter",
                                                         attrs = c(id = paste0("RealParameter.", id$realParameter),
                                                                   spec = "parameter.RealParameter",
                                                                   estimate = "false",
                                                                   lower = "0.0",
                                                                   name = "M",
                                                                   upper = "1.0"),
                                                         value = "1.0")),
                                 xmlNode("parameter", attrs = c(estimate = "false",
                                                                id = paste0("ucldMean.c:", x$partition[i]),
                                                                name = "clock.rate",
                                                                spec = "parameter.RealParameter"),
                                         value = "1.0")))
      id$realParameter <- id$realParameter + 1
      assign("counter", id, envir = x$environment)
    }
    if (x$tree == "Fossilized Birth Death"){
      branchRateModel <- 
        xmlNode("branchRateModel",
                attrs = c(id = paste0("RelaxedClock.c:", x$partition[i]),
                          spec = "beast.evolution.branchratemodel.UCRelaxedClockModel",
                          clock.rate = paste0("@ucldMean.c:", x$partition[i]),
                          rateCategories = paste0("@rateCategories.c:", x$partition[i]),
                          tree = paste0("@Tree.t:", x$partition[i])),
                .children = list(xmlNode("LogNormal",
                                         attrs = c(id = paste0("LogNormalDistributionModel.c:", x$partition[i]),
                                                   S = paste0("@ucldStdev.c:", x$partition[i]),
                                                   meanInRealSpace = "true",
                                                   name = "distr"),
                                         value = xmlNode("parameter",
                                                         attrs = c(id = realParameter(x),
                                                                   spec = "parameter.RealParameter",
                                                                   estimate = "false",
                                                                   lower = "0.0",
                                                                   name = "M",
                                                                   upper = "1.0"),
                                                         value = "1.0"))))
    }
  }
  if (x$clock == "Random Local Clock") {
    if (x$tree %in% c("Yule", "Calibrated Yule", "Birth Death")){
      meanClockRate.c <- xmlNode("parameter",
                                 attrs = c(estimate = "false",
                                           id = paste0("meanClockRate.c:", x$partition[i]),
                                           name = "clock.rate",
                                           spec = "parameter.RealParameter"),
                                 value = "1.0")
      branchRateModel <-
        xmlNode("branchRateModel",
                attrs = c(id = paste0("RandomLocalClock.c:", x$partition[i]),
                          indicators = paste0("@Indicators.c:", x$partition[i]),
                          rates = paste0("@clockrates.c:", x$partition[i]),
                          spec = "beast.evolution.branchratemodel.RandomLocalClockModel",
                          tree = paste0("@Tree.t:", x$partition[i])),
                .children = list(meanClockRate.c))
    } else {
      ## Fossilized Birth Death
      branchRateModel <-
        xmlNode("branchRateModel",
                attrs = c(clock.rate = paste0("@meanClockRate.c:", x$partition[i]),
                          id = paste0("RandomLocalClock.c:", x$partition[i]),
                          indicators = paste0("@Indicators.c:", x$partition[i]),
                          rates = paste0("@clockrates.c:", x$partition[i]),
                          spec = "beast.evolution.branchratemodel.RandomLocalClockModel",
                          tree = paste0("@Tree.t:", x$partition[i])))
    }
  }
  
  # ## branchRateModel
  # ## ---------------
  # if (i == 1){
  #   clocks <- ifelse(link.clocks, "clocks", x$partition[i])
  #   branchRateModel <- xmlNode(
  #     "branchRateModel",
  #     attrs = c(id = paste0("StrictClock.c:", clocks),
  #               spec = "beast.evolution.branchratemodel.StrictClockModel"),
  #     xmlNode("parameter", "1.0", 
  #             attrs = c(estimate = "false",
  #                       id = paste0("clockRate.c:", clocks),
  #                       name = "clock.rate")))
  #   distribution <- addChildren(distribution, kids = list(branchRateModel))
  # } else {
  #   if (!link.clocks){
  #     branchRateModel <- xmlNode(
  #       "branchRateModel",
  #       attrs = c(clock.rate = paste0("@clockRate.c:", x$partition[i]),
  #                 id = paste0("StrictClock.c:", x$partition[i]),
  #                 spec = "beast.evolution.branchratemodel.StrictClockModel"))
  #     distribution <- addChildren(distribution, kids = list(branchRateModel))
  #   }
  # }
  treeLikelihood <- addChildren(treeLikelihood, kids = list(siteModel, branchRateModel))
  likelihood <- addChildren(likelihood, kids = list(treeLikelihood))
  likelihood
}

