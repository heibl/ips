#' @export

operator <- function(ops, x){
  
  i <- 1
  
  ## ID differs for linked or unlinked trees!
  ## ----------------------------------------
  # id <- ifelse(x$link.trees, "trees", x$id[i])
  id <- x$partition[i]
  
  ## Library of possible BEAST parameters (not exhaustive!)
  ## ------------------------------------------------------
  z <- list(
    
    ## Relaxed Clock Exponential
    ## -------------------------
    
    ExpCategoriesRandomWalk.c = list(attrs = c(id = paste0("ExpCategoriesRandomWalk.c:", id),
                                               parameter = paste0("@expRateCategories.c:", id),
                                               spec = "IntRandomWalkOperator",
                                               weight = "10.0",
                                               windowSize = "1")),
    ExpCategoriesSwapOperator.c = list(attrs = c(id = paste0("ExpCategoriesSwapOperator.c:", id),
                                                 intparameter = paste0("@expRateCategories.c:", id),
                                                 spec = "SwapOperator",
                                                 weight = "10.0")),
    ExpCategoriesUniform.c = list(attrs = c(id = paste0("ExpCategoriesUniform.c:", id),
                                            parameter = paste0("@expRateCategories.c:", id),
                                            spec = "UniformOperator",
                                            weight = "10.0")),
    
    ## Relaxed Clock Exponential (and Fossilized Birth Death)
    ## ------------------------------------------------------
    ucedMeanScaler.c = list(attrs = c(id = paste0("ucedMeanScaler.c:", id),
                                      spec = "ScaleOperator",
                                      parameter = paste0("@ucedMean.c:", id),
                                      scaleFactor = "0.5",
                                      weight = "1.0")),
    relaxedUpDownOperatorExp.c = list(attrs = c(id = paste0("relaxedUpDownOperatorExp.c:", id),
                                                spec = "UpDownOperator",
                                                scaleFactor = "0.75",
                                                weight = "3.0"),
                                      value = list(xmlNode("up", attrs = c(idref = paste0("ucedMean.c:", id))),
                                                   xmlNode("down", attrs = c(idref = paste0("Tree.t:", id))))),
    
    ## Relaxed Clock Log Normal
    ## -------------------------
    
    ucldMeanScaler.c = list(attrs = c(id = paste0("ucldMeanScaler.c:", id),
                                      parameter = paste0("@ucldMean.c:", id),
                                      scaleFactor = "0.5",
                                      spec = "ScaleOperator",
                                      weight = "1.0")),
    ucldStdevScaler.c = list(attrs = c(id = paste0("ucldStdevScaler.c:", id),
                                       parameter = paste0("@ucldStdev.c:", id),
                                       scaleFactor = "0.5",
                                       spec = "ScaleOperator",
                                       weight = "3.0")),
    CategoriesRandomWalk.c = list(attrs = c(id = paste0("CategoriesRandomWalk.c:", id),
                                            parameter = paste0("@rateCategories.c:", id),
                                            spec = "IntRandomWalkOperator",
                                            weight = "10.0",
                                            windowSize = "1")),
    CategoriesSwapOperator.c = list(attrs = c(id = paste0("CategoriesSwapOperator.c:", id),
                                              intparameter = paste0("@rateCategories.c:", id),
                                              spec = "SwapOperator",
                                              weight = "10.0")),
    CategoriesUniform.c = list(attrs = c(id = paste0("CategoriesUniform.c:", id),
                                         parameter = paste0("@rateCategories.c:", id),
                                         spec = "UniformOperator",
                                         weight = "10.0")),
    
    relaxedUpDownOperator.c = list(attrs = c(id = paste0("relaxedUpDownOperator.c:", id),
                                             scaleFactor = "0.75",
                                             spec = "UpDownOperator",
                                             weight = "3.0"),
                                   value = list(xmlNode("up", attrs = c(idref = paste0("ucldMean.c:", id))),
                                                xmlNode("down", attrs = c(idref = paste0("Tree.t:", id))))),
    
    ## Random Local Clock
    ## -------------------------
    randomClockScaler.c = list(attrs = c(id = paste0("randomClockScaler.c:", id),
                                      parameter = paste0("@meanClockRate.c:", id),
                                      scaleFactor = "0.5",
                                      spec = "ScaleOperator",
                                      weight = "1.0")),
    randomClockUpDownOperator.c = list(attrs = c(id = paste0("randomClockUpDownOperator.c:", id),
                                             scaleFactor = "0.75",
                                             spec = "UpDownOperator",
                                             weight = "3.0"),
                                   value = list(xmlNode("up", attrs = c(idref = paste0("meanClockRate.c:", id))),
                                                xmlNode("down", attrs = c(idref = paste0("Tree.t:", id))))),
    IndicatorsBitFlip.c = list(attrs = c(id = paste0("IndicatorsBitFlip.c:", id),
                                         parameter = paste0("@Indicators.c:", id),
                                         spec = "BitFlipOperator",
                                         weight = "15.0")),
    ClockRateScaler.c = list(attrs = c(id = paste0("ClockRateScaler.c:", id),
                                       parameter = paste0("@clockrates.c:", id),
                                       scaleFactor = "0.5",
                                       spec = "ScaleOperator",
                                       weight = "15.0")),
    
    ## Yule
    ## ----
    
    YuleBirthRateScaler.t = list(attrs = c(id = paste0("YuleBirthRateScaler.t:", id),
                                           parameter = paste0("@birthRate.t:", id),
                                           scaleFactor = "0.75",
                                           spec = "ScaleOperator",
                                           weight = "3.0")),
    YuleModelTreeScaler.t = list(attrs = c(id = paste0("YuleModelTreeScaler.t:", id),
                                           scaleFactor = "0.5",
                                           spec = "ScaleOperator",
                                           tree = paste0("@Tree.t:", id),
                                           weight = "3.0")),
    YuleModelTreeRootScaler.t = list(attrs = c(id = paste0("YuleModelTreeRootScaler.t:", id),
                                               rootOnly = "true",
                                               scaleFactor = "0.5",
                                               spec = "ScaleOperator",
                                               tree = paste0("@Tree.t:", id),
                                               weight = "3.0")),
    YuleModelUniformOperator.t = list(attrs = c(id = paste0("YuleModelUniformOperator.t:", id),
                                                spec = "Uniform",
                                                tree = paste0("@Tree.t:", id),
                                                weight = "30.0")),
    YuleModelSubtreeSlide.t = list(attrs = c(id = paste0("YuleModelSubtreeSlide.t:", id),
                                             spec = "SubtreeSlide",
                                             tree = paste0("@Tree.t:", id),
                                             weight = "15.0")),
    YuleModelNarrow.t = list(attrs = c(id = paste0("YuleModelNarrow.t:", id),
                                       spec = "Exchange",
                                       tree = paste0("@Tree.t:", id),
                                       weight = "15.0")),
    YuleModelWide.t = list(attrs = c(id = paste0("YuleModelWide.t:", id),
                                     isNarrow = "false",
                                     spec = "Exchange",
                                     tree = paste0("@Tree.t:", id),
                                     weight = "3.0")),
    YuleModelWilsonBalding.t = list(attrs = c(id = paste0("YuleModelWilsonBalding.t:", id),
                                              spec = "WilsonBalding",
                                              tree = paste0("@Tree.t:", id),
                                              weight = "3.0")),
    
    ## Birth Death
    ## -----------
    BirthDeathTreeScaler.t = list(attrs = c(id = paste0("BirthDeathTreeScaler.t:", id),
                                            spec = "ScaleOperator",
                                            scaleFactor = "0.5",
                                            tree = paste0("@Tree.t:", id),
                                            weight = "3.0")),
    BirthDeathTreeRootScaler.t = list(attrs = c(id = paste0("BirthDeathTreeRootScaler.t:", id),
                                                spec = "ScaleOperator",
                                                rootOnly = "true",
                                                scaleFactor = "0.5",
                                                tree = paste0("@Tree.t:", id),
                                                weight = "3.0")),
    BirthDeathUniformOperator.t = list(attrs = c(id = paste0("BirthDeathUniformOperator.t:", id),
                                                 spec = "Uniform",
                                                 tree = paste0("@Tree.t:", id),
                                                 weight = "30.0")),
    BirthDeathSubtreeSlide.t = list(attrs = c(id = paste0("BirthDeathSubtreeSlide.t:", id),
                                              spec = "SubtreeSlide",
                                              tree = paste0("@Tree.t:", id),
                                              weight = "15.0")),
    BirthDeathNarrow.t = list(attrs = c(id = paste0("BirthDeathNarrow.t:", id),
                                        spec = "Exchange",
                                        tree = paste0("@Tree.t:", id),
                                        weight = "15.0")),
    BirthDeathWide.t = list(attrs = c(id = paste0("BirthDeathWide.t:", id),
                                      spec = "Exchange",
                                      isNarrow = "false",
                                      tree = paste0("@Tree.t:", id),
                                      weight = "3.0")),
    BirthDeathWilsonBalding.t = list(attrs = c(id = paste0("BirthDeathWilsonBalding.t:", id),
                                               spec = "WilsonBalding",
                                               tree = paste0("@Tree.t:", id),
                                               weight = "3.0")),
    BirthRateScaler.t = list(attrs = c(id = paste0("BirthRateScaler.t:", id),
                                       spec = "ScaleOperator",
                                       parameter = paste0("@BDBirthRate.t:", id),
                                       scaleFactor = "0.75",
                                       weight = "3.0")),
    DeathRateScaler.t = list(attrs = c(id = paste0("DeathRateScaler.t:", id),
                                       spec = "ScaleOperator",
                                       parameter = paste0("@BDDeathRate.t:", id),
                                       scaleFactor = "0.75",
                                       weight = "3.0")),
    
    ## Strict Clock & Fossilized Birth Death
    ## -------------------------------------
    
    originScalerFBD.t = list(attrs = c(id = paste0("originScalerFBD.t:", id),
                                       parameter = paste0("@originFBD.t:", id),
                                       scaleFactor = "0.75",
                                       spec = "ScaleOperator",
                                       weight = "3.0")),
    divRateScalerFBD.t = list(attrs = c(id = paste0("divRateScalerFBD.t:", id),
                                        parameter = paste0("@diversificationRateFBD.t:", id),
                                        scaleFactor = "0.75",
                                        spec = "ScaleOperator",
                                        weight = "10.0")),
    turnoverScalerFBD.t = list(attrs = c(id = paste0("turnoverScalerFBD.t:", id),
                                         parameter = paste0("@turnoverFBD.t:", id),
                                         scaleFactor = "0.75",
                                         spec = "ScaleOperator",
                                         weight = "10.0")),
    samplingPScalerFBD.t = list(attrs = c(id = paste0("samplingPScalerFBD.t:", id),
                                          parameter = paste0("@samplingProportionFBD.t:", id),
                                          scaleFactor = "0.75",
                                          spec = "ScaleOperator",
                                          weight = "10.0")),
    LeafToSAFBD.t = list(attrs = c(id = paste0("LeafToSAFBD.t:", id),                                   
                                   spec = "LeafToSampledAncestorJump",
                                   tree = paste0("@Tree.t:", id),
                                   weight = "10.0")),
    SAWilsonBaldingFBD.t = list(attrs = c(id = paste0("SAWilsonBaldingFBD.t:", id), 
                                          spec = "SAWilsonBalding",
                                          tree = paste0("@Tree.t:", id),
                                          weight = "10.0")),
    SAWideFBD.t = list(attrs = c(id = paste0("SAWideFBD.t:", id),
                                 isNarrow = "false",
                                 spec = "SAExchange",
                                 tree = paste0("@Tree.t:", id),
                                 weight = "10.0")),
    SANarrowFBD.t = list(attrs = c(id = paste0("SANarrowFBD.t:", id),
                                   spec = "SAExchange",
                                   tree = paste0("@Tree.t:", id),
                                   weight = "10.0")),
    SAUniformOperatorFBD.t = list(attrs = c(id = paste0("SAUniformOperatorFBD.t:", id),
                                            spec = "SAUniform",
                                            tree = paste0("@Tree.t:", id),
                                            weight = "20.0")),
    SATreeRootScalerFBD.t = list(attrs = c(id = paste0("SATreeRootScalerFBD.t:", id),
                                           rootOnly = "true",
                                           scaleFactor = "0.95",
                                           spec = "SAScaleOperator",
                                           tree = paste0("@Tree.t:", id),
                                           weight = "1.0")),
    SATreeScalerFBD.t = list(attrs = c(id = paste0("SATreeScalerFBD.t:", id),
                                       # rootOnly = "true", ## not in UCL + FBD
                                       scaleFactor = "0.95",
                                       spec = "SAScaleOperator",
                                       tree = paste0("@Tree.t:", id),
                                       weight = "3.0")),
    StrictClockRateScaler.c = list(attrs = c(id = paste0("StrictClockRateScaler.c:", id),
                                             parameter = paste0("@clockRate.c:", id),
                                             scaleFactor = "0.75",
                                             spec = "ScaleOperator",
                                             weight = "3.0")),
    strictClockUpDownOperator.c = list(attrs = c(id = paste0("strictClockUpDownOperator.c:", id),
                                                 scaleFactor = "0.75",
                                                 spec = "UpDownOperator",
                                                 weight = "3.0"),
                                       value = list(xmlNode("up", attrs = c(idref = paste0("clockRate.c:", id))),
                                                    xmlNode("down", attrs = c(idref = paste0("Tree.t:", id))))),
    ## Substitution models
    
    KappaScaler.s = list(attrs = c(id = paste0("KappaScaler.s:", id),
                                   parameter = paste0("@kappa.s:", id),
                                   scaleFactor = "0.5",
                                   spec = "ScaleOperator",
                                   weight = "0.1")),
    kappa1Scaler.s = list(attrs = c(id = paste0("kappa1Scaler.s:", id),
                                    parameter = paste0("@kappa1.s:", id),
                                    scaleFactor = "0.5",
                                    spec = "ScaleOperator",
                                    weight = "0.1")),
    kappa2Scaler.s = list(attrs = c(id = paste0("kappa2Scaler.s:", id),
                                    parameter = paste0("@kappa2.s:", id),
                                    scaleFactor = "0.5",
                                    spec = "ScaleOperator",
                                    weight = "0.1")),
    RateACScaler.s = list(attrs = c(id = paste0("RateACScaler.s:", id),
                                    parameter = paste0("@rateAC.s:", id),
                                    scaleFactor = "0.5",
                                    spec = "ScaleOperator",
                                    weight = "0.1")),
    RateAGScaler.s = list(attrs = c(id = paste0("RateAGScaler.s:", id),
                                    parameter = paste0("@rateAG.s:", id),
                                    scaleFactor = "0.5",
                                    spec = "ScaleOperator",
                                    weight = "0.1")),
    RateATScaler.s = list(attrs = c(id = paste0("RateATScaler.s:", id),
                                    parameter = paste0("@rateAT.s:", id),
                                    scaleFactor = "0.5",
                                    spec = "ScaleOperator",
                                    weight = "0.1")),
    RateCGScaler.s = list(attrs = c(id = paste0("RateCGScaler.s:", id),
                                    parameter = paste0("@rateCG.s:", id),
                                    scaleFactor = "0.5",
                                    spec = "ScaleOperator",
                                    weight = "0.1")),
    RateGTScaler.s = list(attrs = c(id = paste0("RateGTScaler.s:", id),
                                    parameter = paste0("@rateGT.s:", id),
                                    scaleFactor = "0.5",
                                    spec = "ScaleOperator",
                                    weight = "0.1")),
    FrequenciesExchanger.s = list(attrs = c(delta = "0.01",
                                            id = paste0("FrequenciesExchanger.s:", id),
                                            spec = "DeltaExchangeOperator",
                                            weight = "0.1"),
                                  value = list(xmlNode("parameter", 
                                                       attrs = c(idref = paste0("freqParameter.s:", id)))))
    
    #   +     <operator delta='0.01' id='FrequenciesExchanger.s:bears' spec='DeltaExchangeOperator' weight='0.1'>
    #   +       <parameter idref='freqParameter.s:bears' />
    #   +     </operator>
  )
  
  
  
  ops <- match.arg(ops, names(z))
  
  z <- z[[ops]]
  if (is.null(z$value)){
    xmlNode("operator", attrs = z$attrs)
  } else {
    xmlNode("operator", attrs = z$attrs, .children = z$value)
  }
}

# operator_list <- read.csv("data/operator_list.csv", stringsAsFactors = FALSE)
# save(operator_list, file = "data/operator_list.rda")