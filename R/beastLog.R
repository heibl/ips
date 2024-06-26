## This code is part of the ips package
## © C. Heibl 2019 (last update 2020-03-13)

#' @title XML Parameter Nodes
#' @description A collection of parameter nodes used in the standard template of BEAST2.
#' @param l A character string giving the name (i.e. \code{id} or \code{idref} attribute) of the desired parameter.
#' @param i Integer giving the index to a partition, internally generated by \code{rbeauti}.
#' @param x A list, internally generated by \code{rbeauti}.
#' @return An object of class \code{XMLNode}.
#' @export
#' @keywords internal

beastLog <- function(l, x, i){
  
  ## ID differs for linked or unlinked trees!
  ## ----------------------------------------
  # id <- ifelse(x$link.trees, "trees", x$partition[i])
  id <- x$partition[i]
  
  ## Library of possible BEAST parameters (not exhaustive!)
  ## ------------------------------------------------------
  z <- list(
    
    posterior = list(attrs = c(idref = "posterior")),
    likelihood = list(attrs = c(idref = "likelihood")),
    prior = list(attrs = c(idref = "prior")),
    treeLikelihood = list(attrs = c(idref = paste0("treeLikelihood.", id))),
    TreeHeight.t = list(attrs = c(id = paste0("TreeHeight.t:", id),
                                  spec = "beast.evolution.tree.TreeHeightLogger",
                                  tree = paste0("@Tree.t:", id))),
    
    
    
    
    ## Strict Clock
    clockRate.c = list(attrs = c(idref = paste0("clockRate.c:", id))), 
    
    ## Relaxed Clock Exponential
    rateStat.c = list(attrs = c(branchratemodel = paste0("@ExponentialRelaxedClock.c:", id),
                                id = paste0("rateStat.c:", id),
                                spec = "beast.evolution.branchratemodel.RateStatistic",
                                tree = paste0("@Tree.t:", id))),
    ucedMean.c = list(attrs = c(idref = paste0("ucedMean.c:", id))),
    
    ## Relaxed Clock Log Normal
    ucldMean.c = list(attrs = c(idref = paste0("ucldMean.c:", id))),
    ucldStdev.c = list(attrs = c(idref = paste0("ucldStdev.c:", id))),
    rate.c = list(attrs = c(branchratemodel = paste0("@RelaxedClock.c:", id),
                            id = paste0("rate.c:", id),
                            spec = "beast.evolution.branchratemodel.RateStatistic",
                            tree = paste0("@Tree.t:", id))),
    
    ## Random Local Clock
    Indicators.c = list(attrs = c(idref = paste0("Indicators.c:", id))),
    clockrates.c = list(attrs = c(idref = paste0("clockrates.c:", id))),
    meanClockRate.c = list(attrs = c(idref = paste0("meanClockRate.c:", id))),
    RRateChanges.c = list(attrs = c(idref = paste0("RRateChanges.c:", id))),
    
    ## Yule
    YuleModel.t = list(attrs = c(idref = paste0("YuleModel.t:", id))),
    birthRate.t = list(attrs = c(idref = paste0("birthRate.t:", id))),
    
    ## Birth Death
    BirthDeath.t = list(attrs = c(idref = paste0("BirthDeath.t:", id))),
    BDBirthRate.t = list(attrs = c(idref = paste0("BDBirthRate.t:", id))),
    BDDeathRate.t = list(attrs = c(idref = paste0("BDDeathRate.t:", id))),
    
    ## Fossilized Birth Death
    FBD.t = list(attrs = c(idref = paste0("FBD.t:", id))),
    diversificationRateFBD.t = list(attrs = c(idref = paste0("diversificationRateFBD.t:", id))),
    turnoverFBD.t = list(attrs = c(idref = paste0("turnoverFBD.t:", id))),
    samplingProportionFBD.t = list(attrs = c(idref = paste0("samplingProportionFBD.t:", id))),
    originFBD.t = list(attrs = c(idref = paste0("originFBD.t:", id))),
    SACountFBD.t = list(attrs = c(id = paste0("SACountFBD.t:", id),
                                  spec = "beast.evolution.tree.SampledAncestorLogger",
                                  tree = paste0("@Tree.t:", id))),
    
    ## HKY
    kappa.s = list(attrs = c(idref = paste0("kappa.s:", id))),
    kappa1.s = list(attrs = c(idref = paste0("kappa1.s:", id))),
    kappa2.s = list(attrs = c(idref = paste0("kappa2.s:", id))),
    rateAC.s = list(attrs = c(idref = paste0("rateAC.s:", id))),
    rateAG.s = list(attrs = c(idref = paste0("rateAG.s:", id))),
    rateAT.s = list(attrs = c(idref = paste0("rateAT.s:", id))),
    rateCG.s = list(attrs = c(idref = paste0("rateCG.s:", id))),
    rateGT.s = list(attrs = c(idref = paste0("rateGT.s:", id))),
    freqParameter.s = list(attrs = c(idref = paste0("freqParameter.s:", id)))
    
  )
  l <- match.arg(l, names(z))
  
  z <- z[[l]]
  if (is.null(z$value)){
    xmlNode("log", attrs = z$attrs)
  } else {
    xmlNode("log", attrs = z$attrs, .children = z$value)
  }
}

# log_list <- read.csv2("data/log_list.csv", stringsAsFactors = FALSE)
# save(log_list, file = "data/log_list.rda")