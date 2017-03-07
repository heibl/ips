## This code is part of the ips package
## © C. Heibl 2014 (last update 2015-04-03)

setMCMC  <- function(xml, chainLength, storeEvery, pre){
  
  if ( !missing (chainLength) ){
    n <- xmlRoot(xml)[["run"]]
    xmlAttrs(n) <- c(chainLength = chainLength)
    replaceNodes(xmlRoot(xml)[["run"]], n)
  }
  if ( !missing (storeEvery) ){
    n <- xmlRoot(xml)[["run"]][["state"]]
    xmlAttrs(n) <- c(storeEvery = storeEvery)
    replaceNodes(xmlRoot(xml)[["run"]][["state"]], n)
  }
  xml
}
