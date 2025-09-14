## This code is part of the ips package
## © C. Heibl 2014 (last update 2025-08-30)

#' @importFrom XML addChildren

assembleLoggers <- function(id){
  
  ## assemble node <logger id="tracelog">
  ## ------------------------------------
  tracelog <- xmlNode("logger", 
                      attrs = c(fileName = paste(id[1], ".log", sep = ""),
                                          id = "tracelog",
                                          logEvery = "1000",
                                          model = "@posterior",
                                          sanitiseHeaders = "true",
                                          sort="smart"),
                      .children = list(
                        xmlNode("log", attrs = c(idref = "posterior")),
                        xmlNode("log", attrs = c(idref = "likelihood")),
                        xmlNode("log", attrs = c(idref = "prior"))))
  
  ## assemble node <logger id="screenlog">
  ## -------------------------------------
  screenlog <- xmlNode("logger", 
                       attrs = c(id = "screenlog", logEvery = "1000"),
                       .children = list(xmlNode("log", attrs = c(idref = "posterior")),
                         xmlNode("log", attrs = c(arg = "@posterior",
                                                  id = "ESS.0",
                                                  spec = "util.ESS")),
                         xmlNode("log", attrs = c(idref = "likelihood")),
                         xmlNode("log", attrs = c(idref = "prior"))))
  
  treelog <- vector(mode = "list")
  
  for ( i in seq_along(id) ){
    
    ## add children to node <logger id="tracelog">
    ## -------------------------------------------
    tracelog <- addChildren(tracelog, kids = list(
      xmlNode("log", attrs = c(idref = paste("treeLikelihood.", id[i], sep = ""))),
      xmlNode("log", attrs = c(id = paste("TreeHeight.t:", id[i], sep = ""),
                               spec = "beast.evolution.tree.TreeHeightLogger",
                               tree = paste("@Tree.t:", id[i], sep = "")))))
    if ( i > 1 ){
      tracelog <- addChildren(tracelog, kids = list(
        xmlNode("parameter", attrs = c(idref = paste("clockRate.c:", id[i], sep = ""), 
                                       name = "log"))))
    }
    tracelog <- addChildren(tracelog, kids = list(
      xmlNode("log", attrs = c(idref = paste("YuleModel.t:", id[i], sep = ""))),
      xmlNode("parameter", attrs = c(idref = paste("birthRate.t:", id[i], sep = ""), 
                                     name = "log"))))
    
#     if ( !missing(taxonset) ){
#       tslog <- lapply(taxonset, function(x) xmlNode("log", attrs = c(idref = paste(x$id, ".prior", sep = ""))))
#       logger1 <- addChildren(logger1, kids = tslog)
#     }

    ## treelog
    ## -------
    this.treelog <- xmlNode("logger", 
                            attrs = c(fileName = "$(tree).trees",
                                                id = paste("treelog.t:", id[i], sep = ""),
                                                logEvery = "1000",
                                                mode = "tree"),
                            .children = list(xmlNode("log", 
                                                     attrs = c(id = paste("TreeWithMetaDataLogger.t:", id[i], sep = ""),
                                                                      spec = "beast.evolution.tree.TreeWithMetaDataLogger",
                                                                      tree = paste("@Tree.t:", id[i], sep = "")))))
    treelog <- c(treelog, list(this.treelog))
  }
  c(list(tracelog), list(screenlog), treelog)
}