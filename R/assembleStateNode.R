## This code is part of the ips package
## © C. Heibl 2014 (last update 2015-04-04)

#' @importFrom XML addChildren

assembleStateNode <- function(id){
  
  state <- xmlNode("state", 
                   attrs = c(id = "state",
                             storeEvery = "5000"))
  for ( i in seq_along(id) ){
    
    ## tree: taxonset
    ## --------------
    taxonset <- xmlNode("taxonset", 
                          attrs = c(id = paste("TaxonSet.", id[i], sep = ""), spec = "TaxonSet"),
                          .children = list(xmlNode("data",
                                                   attrs = c(idref = id[i],
                                                             name = "alignment"))))
    tree <- xmlNode("tree", 
                    attrs = c(id = paste("Tree.t:", id[i], sep = ""),
                              name = "stateNode"),
                    .children = list(taxonset))
    
    state <- addChildren(state, kids = list(tree))
    
    ## parameter: clockRate
    ## --------------------
    if ( i > 1 ){
      br <- xmlNode("parameter", "1.0",
                    attrs = c(id = paste("clockRate.c:", id[i], sep = ""),
                              name = "stateNode"))
      state <- addChildren(state, kids = list(br))
    }
    
    ## parameter: birthRate
    ## --------------------
    br <- xmlNode("parameter", "1.0",
            attrs = c(id = paste("birthRate.t:", id[i], sep = ""),
                      name = "stateNode"))
    state <- addChildren(state, kids = list(br))
    
  }
  state
}