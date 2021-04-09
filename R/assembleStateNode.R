## This code is part of the ips package
## © C. Heibl 2014 (last update 2020-03-13)

## To do: stateNodes must be added

#' @export

assembleStateNode <- function(x){
  
  state <- xmlNode("state", 
                   attrs = c(id = "state",
                             spec = "State",
                             storeEvery = "5000"))
  
  
  ## Loop over partitions
  ## --------------------
  for (i in seq_along(x$partition)){
    
    ## Tree.t
    if (i == 1 | !x$link.trees){
      state <- addChildren(state, kids = list(assembleTreeNode(i, x)))
    }
    

    ## Substitution models  
    ## -------------------
    if (x$subst.model == "HKY"){
      params <- c("kappa.s", "freqParameter.s")
      kids <- lapply(params, parameter, i = i, x = x)
      state <- addChildren(state, kids = kids)
    }
    if (x$subst.model == "TN93"){
      params <- c("kappa1.s", "kappa2.s", "freqParameter.s")
      kids <- lapply(params, parameter, i = i, x = x)
      state <- addChildren(state, kids = kids)
    }
    if (x$subst.model == "GTR"){
      params <- c("rateAC.s", "rateAG.s", "rateAT.s", "rateGC.s", "rateGT.s", "freqParameter.s")
      kids <- lapply(params, parameter, i = i, x = x)
      state <- addChildren(state, kids = kids)
    }
    
    ## RANDOM LOCAL CLOCK  
    ## ------------------
    if (x$clock == "Random Local Clock"){
      
      # if (is.na(x$tip.dates)){
      #   params <- "ucldStdev.c"
      # } else {
      if (x$tree %in% c("Yule", "Calibrated Yule", "Birth Death")){
        params <- c("Indicators.c", "clockrates.c")
      }
      if (x$tree == "Fossilized Birth Death"){
        params <- c("Indicators.c", "meanClockRate.c", "clockrates.c")
      }
      # }
      kids <- lapply(params, parameter, i = i, x = x)
      state <- addChildren(state, kids = kids)
    }
    
    ## *.t 
    ## Choose parameters according to tree prior
    if (x$tree == "Yule"){
      par_set <- "birthRate.t"
    }
    if (x$tree == "Calibrated Yule"){
      par_set <- "birthRateY.t"
    }
    if (x$tree == "Birth Death"){
      par_set <- c("BDBirthRate.t", "BDDeathRate.t")
    }
    if (x$tree == "Fossilized Birth Death"){
      par_set <- c("diversificationRateFBD.t", "turnoverFBD.t", "samplingProportionFBD.t", "originFBD.t") 
    }
    if (i == 1 | !x$link.trees){
      state <- addChildren(state, kids = lapply(par_set, parameter, i = i, x = x))
    }
  
    ## clockRate.c (only add if clock models are unlinked)
    ## ---------------------------------------------------
    # if ((prior != "Yule" & clock != "Strict Clock" ) | (!link.clocks & i > 1)){    ## needs to be updated
    
    ## Choose parameters according to tree prior
    if (x$clock == "Strict Clock"){
      # clock_par <- "clockRate.c" ## should be added when i > 1
      if (x$tree %in% c("Yule", "Calibrated Yule", "Birth Death")) {
        clock_par <- NULL
      } else {
        clock_par = "clockRate.c"
      }
      state <- addChildren(state, kids = lapply(clock_par, parameter, i = i, x = x))
    }
    if (x$clock == "Relaxed Clock Exponential"){
      if (x$tree %in% c("Yule", "Birth Death")) {
        kids <- list(xmlNode("stateNode",
                             attrs = c(dimension = 4,
                                       id = paste0("expRateCategories.c:", x$partition[i]),
                                       spec = "parameter.IntegerParameter"),
                             .children = 1))
        ## ucedMean.c defined at branchratemodel
      } else {
        kids <- list(xmlNode("stateNode",
                             attrs = c(dimension = 4,
                                       id = paste0("expRateCategories.c:", x$partition[i]),
                                       spec = "parameter.IntegerParameter"),
                             .children = 1),
                     parameter("ucedMean.c.FBD", i = i, x = x))
      }
      state <- addChildren(state, kids = kids)
    }
    
    if (x$clock == "Relaxed Clock Log Normal"){
      
      if (is.na(x$tip.dates)){
        params <- "ucldStdev.c"
      } else {
        params <- c("ucldMean.c", "ucldStdev.c")
      }
      kids <- lapply(params, parameter, i = i, x = x)
      kids <- c(kids,
                list(xmlNode("stateNode",
                             attrs = c(dimension = 4,
                                       id = paste0("rateCategories.c:", x$partition[i]),
                                       spec = "parameter.IntegerParameter"),
                             .children = 1)))
      state <- addChildren(state, kids = kids)
    }

    
  
    
    
    
  }
  state
}
