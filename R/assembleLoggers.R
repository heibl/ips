## This code is part of the ips package
## © C. Heibl 2014 (last update 2020-03-32)

#' @importFrom utils data
#' @export

assembleLoggers <- function(x, taxonset){
  
  data(log_list, envir = environment())
  i <- 1
  id <- x$partition[i]
  
  ## Assemble node <logger id="tracelog">
  ## ------------------------------------
  if (x$tree != "Fossilized Birth Death"){
    log_list <- log_list[log_list$tip.dates != is.na(x$tip.dates), ]
  }
  params <- log_list$name[log_list$clock == x$clock & log_list$tree == x$tree]
  if (x$subst.model == "HKY") params <- c("kappa.s", params, "freqParameter.s")
  if (x$subst.model == "TN93") params <- c("kappa1.s", "kappa2.s", params, "freqParameter.s")
  if (x$subst.model == "GTR") params <- c("rateAC.s", "rateAG.s", "rateAT.s", 
                                          "rateCG.s", "rateGT.s", params, "freqParameter.s")
  params <- c("posterior", "likelihood", "prior", "treeLikelihood", "TreeHeight", params)
  
  tracelog <- xmlNode("logger", 
                      attrs = c(fileName = x$file.name$log,
                                id = "tracelog",
                                logEvery = x$log.every,
                                model = "@posterior",
                                sanitiseHeaders = "true",
                                sort = "smart",
                                spec = "Logger"),
                      .children = lapply(params, beastLog, x = x, i = i))
  
  ## Assemble node <logger id="screenlog">
  ## -------------------------------------
  screenlog <- xmlNode("logger", 
                       attrs = c(id = "screenlog", 
                                 logEvery = x$log.every,
                                 spec = "Logger"),
                       .children = lapply(c("posterior", "likelihood", "prior"),
                                          beastLog, x = x, i = i))
  
  ## Assemble node <logger id="treelog">
  ## -------------------------------------
  treelog <- vector(mode = "list")
  
  
  ## Add *tree*-related children ... 
  ## -------------------------------
  for (i in seq_along(x$partition)){
    
    ## ... for *linked* trees ... 
    ## --------------------------
    if (x$link.trees){
      if (i == 1){
        
        ## ... to node <logger id="tracelog">
        ## ----------------------------------
        # tracelog <- addChildren(tracelog, kids = list(
        #   xmlNode("log", attrs = c(id = paste0("TreeHeight.t:", id),
        #                            spec = "beast.evolution.tree.TreeHeightLogger",
        #                            tree = paste0("@Tree.t:", id))),
        #   xmlNode("log", attrs = c(idref = paste0("YuleModel.t:", id))),
        #   xmlNode("parameter", attrs = c(idref = paste0("birthRate.t:", id),
        #                                  name = "log"))))
        
        if (!missing(taxonset)){
          tslog <- lapply(names(taxonset), function(x) xmlNode("log", 
                                                               attrs = c(idref = paste(x, "trees.prior", sep = "."))))
          tracelog <- addChildren(tracelog, kids = tslog)
        }
        
        ## ... to node<logger id="treelog">
        ## --------------------------------
        if (x$clock == "Strict Clock"){
          this.treelog <- xmlNode("log", 
                                  attrs = c(id = paste0("TreeWithMetaDataLogger.t:", id),
                                            spec = "beast.evolution.tree.TreeWithMetaDataLogger",
                                            tree = paste0("@Tree.t:", id)))
        }
        if (x$clock == "Relaxed Clock Exponential"){
          this.treelog <- xmlNode("log", 
                                  attrs = c(branchratemodel = paste0("@ExponentialRelaxedClock.c:", id),
                                            id = paste0("TreeWithMetaDataLogger.t:", id),
                                            spec = "beast.evolution.tree.TreeWithMetaDataLogger",
                                            tree = paste0("@Tree.t:", id)))
        }
        if (x$clock == "Relaxed Clock Log Normal"){
          this.treelog <- xmlNode("log", 
                                  attrs = c(branchratemodel = paste0("@RelaxedClock.c:", id),
                                            id = paste0("TreeWithMetaDataLogger.t:", id),
                                            spec = "beast.evolution.tree.TreeWithMetaDataLogger",
                                            tree = paste0("@Tree.t:", id)))
        }
        if (x$clock == "Random Local Clock"){
          this.treelog <- xmlNode("log", 
                                  attrs = c(branchratemodel = paste0("@RandomLocalClock.c:", id),
                                            id = paste0("TreeWithMetaDataLogger.t:", id),
                                            spec = "beast.evolution.tree.TreeWithMetaDataLogger",
                                            tree = paste0("@Tree.t:", id)))
        }
        this.treelog <- xmlNode("logger", 
                                attrs = c(fileName = x$file.name$trees,
                                          id = paste0("treelog.t:", id),
                                          logEvery = x$log.every,
                                          mode = "tree",
                                          spec = "Logger"),
                                .children = list(this.treelog))
        treelog <- c(treelog, list(this.treelog))
        
      } else {
        ## Do nothing
      }
    } else {
      
      ## ... for *unlinked* trees to node <logger id="tracelog"> 
      ## -------------------------------------------------------
      tracelog <- addChildren(tracelog, kids = list(
        xmlNode("log", attrs = c(id = paste0("TreeHeight.t:", x$partition[i]),
                                 spec = "beast.evolution.tree.TreeHeightLogger",
                                 tree = paste0("@Tree.t:", x$partition[i]))),
        xmlNode("log", attrs = c(idref = paste0("YuleModel.t:", x$partition[i]))),
        xmlNode("parameter", attrs = c(idref = paste0("birthRate.t:", x$partition[i]), 
                                       name = "log"))))
      
      if (!missing(taxonset)){
        tslog <- lapply(names(taxonset), function(x) xmlNode("log", 
                                                             attrs = c(idref = paste(x, x$partition[i], "prior", sep = "."))))
        tracelog <- addChildren(tracelog, kids = tslog)
      }
      
      ## ... for *unlinked* trees to node <logger id="treelog"> 
      ## -------------------------------------------------------
      this.treelog <- xmlNode("log", 
                              attrs = c(id = paste0("TreeWithMetaDataLogger.t:", x$partition[i]),
                                        spec = "beast.evolution.tree.TreeWithMetaDataLogger",
                                        tree = paste0("@Tree.t:", x$partition[i])))
      this.treelog <- xmlNode("logger", 
                              attrs = c(fileName = x$file.name$trees,
                                        id = paste0("treelog.t:", x$partition[i]),
                                        logEvery = x$log.every,
                                        mode = "tree",
                                        spec = "Logger"),
                              .children = list(this.treelog))
      treelog <- c(treelog, list(this.treelog))
    }
    
    if (!x$link.clocks & i > 1){
      tracelog <- addChildren(tracelog, kids = list(
        xmlNode("parameter", attrs = c(idref = paste0("clockRate.c:", x$partition[i]), 
                                       name = "log"))))
    }
  }
  c(list(tracelog), list(screenlog), treelog)
}