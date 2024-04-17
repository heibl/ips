## This code is part of the ips package
## © C. Heibl 2020 (last update 2020-02-27)

#' @export

parameter <- function(param, i, x){
  
  ## ID differs for linked or unlinked trees!
  ## ----------------------------------------
  # id <- ifelse(x$link.trees, "trees", x$id[i])
  id <- x$partition[i]
  
  ## Library of possible BEAST parameters (not exhaustive!)
  ## ------------------------------------------------------
  z <- list(
    
    ## Site Model
    mutationRate.s = list(attrs = c(id = paste0("mutationRate.s:", id),
                                    name = "mutationRate",
                                    spec = "parameter.RealParameter",
                                    estimate = "false"),
                          value = "1.0"),
    gammaShape.s = list(attrs = c(id = paste0("gammaShape.s:", id),
                                  name = "shape",
                                  spec = "parameter.RealParameter",
                                  estimate = "false"),
                        value = "1.0"),
    proportionInvariant.s = list(attrs = c(id = paste0("proportionInvariant.s:", id),
                                           name = "proportionInvariant",
                                           spec = "parameter.RealParameter",
                                           lower = "0.0",
                                           upper = "1.0",
                                           estimate = "false"),
                                 value = "0.0"),
    
    ## HKY
    kappa.s = list(attrs = c(id = paste0("kappa.s:", id),
                             lower = "0.0",
                             name = "stateNode",
                             spec = "parameter.RealParameter"),
                   value = "2.0"),
    ## TN93
    kappa1.s = list(attrs = c(id = paste0("kappa1.s:", id),
                              lower = "0.0",
                              name = "stateNode",
                              spec = "parameter.RealParameter"),
                    value = "2.0"),
    kappa2.s = list(attrs = c(id = paste0("kappa2.s:", id),
                              lower = "0.0",
                              name = "stateNode",
                              spec = "parameter.RealParameter"),
                    value = "2.0"),
    ## HKY + TN93
    freqParameter.s = list(attrs = c(dimension = "4",
                                     id = paste0("freqParameter.s:", id),
                                     lower = "0.0",
                                     name = "stateNode",
                                     spec = "parameter.RealParameter",
                                     upper = "1.0"),
                           value = "0.25"),
    ## GTR
    rateAC.s = list(attrs = c(id = paste0("rateAC.s:", id),
                              lower = "0.0",
                              name = "stateNode",
                              spec = "parameter.RealParameter"),
                    value = "1.0"),
    rateAG.s = list(attrs = c(id = paste0("rateAG.s:", id),
                              lower = "0.0",
                              name = "stateNode",
                              spec = "parameter.RealParameter"),
                    value = "1.0"),
    rateAT.s = list(attrs = c(id = paste0("rateAT.s:", id),
                              lower = "0.0",
                              name = "stateNode",
                              spec = "parameter.RealParameter"),
                    value = "1.0"),
    rateGC.s = list(attrs = c(id = paste0("rateCG.s:", id),
                              lower = "0.0",
                              name = "stateNode",
                              spec = "parameter.RealParameter"),
                    value = "1.0"),
    rateGT.s = list(attrs = c(id = paste0("rateGT.s:", id),
                              lower = "0.0",
                              name = "stateNode",
                              spec = "parameter.RealParameter"),
                    value = "1.0"),
    
    clockRate.c = list(attrs = c(id = paste0("clockRate.c:", id),
                                 name = "stateNode",
                                 spec = "parameter.RealParameter"),
                       value = "1.0"),
    UCExpLambda.c = list(attrs = c(id = paste0("UCExpLambda.c:", id),
                                   name = "mean",
                                   spec = "parameter.RealParameter"),
                         value = "1.0"),
    ## Relaxed Clock Exponential / Yule
    ucedMean.c.Yule = list(attrs = c(estimate = "false",
                                     id = paste0("ucedMean.c:", id),
                                     name = "clock.rate",
                                     spec = "parameter.RealParameter"),
                           value = "1.0"),
    ## Relaxed Clock Exponential / Fossilized Birth Death == state node
    ucedMean.c.FBD = list(attrs = c(id = paste0("ucedMean.c:", id),
                                    spec = "parameter.RealParameter",
                                    name = "stateNode"),
                          value = "1.0"),
    ## Relaxed Clock Log Normal
    ucldMean.c = list(attrs = c(id = paste0("ucldMean.c:", id),
                                name = "stateNode",
                                spec = "parameter.RealParameter"),
                      value = "1.0"),
    ucldStdev.c = list(attrs = c(id = paste0("ucldStdev.c:", id),
                                 name = "stateNode",
                                 spec = "parameter.RealParameter",
                                 lower = "0.0"),
                       value = "0.1"),
    
    ## Random Local Clock
    ## ------------------
    Indicators.c = list(attrs = c(dimension = "4",
                                  id = paste0("Indicators.c:", id),
                                  name = "stateNode",
                                  spec = "parameter.BooleanParameter"),
                        value = "false"),
    meanClockRate.c = list(attrs = c(id = paste0("meanClockRate.c:", id),
                                     name = "stateNode",
                                     spec = "parameter.RealParameter"),
                        value = "1.0"),
    clockrates.c = list(attrs = c(dimension = "4",
                                  id = paste0("clockrates.c:", id),
                                  lower = "1.0E-9",
                                  name = "stateNode",
                                  spec = "parameter.RealParameter"),
                        value = "1.0"),
    ## brachRateModel
    # meanClockRate.c = list(attrs = c(estimate = "false",
    #                                  id = paste0("meanClockRate.c:", id),
    #                                  name = "clock.rate",
    #                                  spec = "parameter.RealParameter"),
    #                        value = "1.0"),
    
    ## Yule
    birthRate.t = list(attrs = c(id = paste0("birthRate.t:", id),
                                 name = "stateNode",
                                 spec = "parameter.RealParameter"),
                       value = "1.0"),
    ## Calibrated Yule
    birthRateY.t = list(attrs = c(id = paste0("birthRateY.t:", id),
                                  name = "stateNode",
                                  spec = "parameter.RealParameter"),
                        value = "1.0"),
    ## Birth Death
    BDBirthRate.t = list(attrs = c(id = paste0("BDBirthRate.t:", id),
                                   lower = "0.0",
                                   name = "stateNode",
                                   spec = "parameter.RealParameter",
                                   upper = "10000.0"),
                         value = "1.0"),
    BDDeathRate.t = list(attrs = c(id = paste0("BDDeathRate.t:", id),
                                   lower = "0.0",
                                   name = "stateNode",
                                   spec = "parameter.RealParameter",
                                   upper = "1.0"),
                         value = "0.5"),
    ## Fossilzed Birth Death
    diversificationRateFBD.t = list(attrs = c(id = paste0("diversificationRateFBD.t:", id),
                                              name = "stateNode",
                                              spec = "parameter.RealParameter",
                                              lower = "0.0"),
                                    value = "1.0"),
    turnoverFBD.t = list(attrs = c(id = paste0("turnoverFBD.t:", id),
                                   name = "stateNode",
                                   spec = "parameter.RealParameter",
                                   lower = "0.0",
                                   upper = "1.0"),
                         value = "0.5"),
    samplingProportionFBD.t = list(attrs = c(id = paste0("samplingProportionFBD.t:", id),
                                             name = "stateNode",
                                             spec = "parameter.RealParameter",
                                             lower = "0.0",
                                             upper = "1.0"),
                                   value = "0.5"),
    originFBD.t = list(attrs = c(id = paste0("originFBD.t:", id),
                                 name = "stateNode",
                                 spec = "parameter.RealParameter",
                                 lower = "0.0"),
                       value = "100.0"),
    rFBD.t = list(attrs = c(id = paste0("rFBD.t:", id),
                            name = "removalProbability",
                            spec = "parameter.RealParameter",
                            lower = "0.0",
                            upper = "1.0"),
                  value = "0.0"),
    rhoFBD.t = list(attrs = c(id = paste0("rhoFBD.t:", id),
                              name = "rho",
                              spec = "parameter.RealParameter",
                              estimate = "false",
                              lower = "0.0",
                              upper = "1.0"),
                    value = "1.0")
  )
  
  param <- match.arg(param, names(z))
  
  z <- z[[param]]
  xmlNode("parameter", attrs = z$attrs, .children = list(z$value))
}