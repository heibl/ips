rbeauti.taxonset <- function(x, id){
  
  taxon <- lapply(x$taxon, function(x) xmlNode("taxon", 
                                               attrs = c(id = x,
                                                         spec = "Taxon")))
  taxonset <- xmlNode("taxonset",
                      attrs = c(id = x$id,
                                spec = "TaxonSet"),
                      .children = taxon)
  
  LogNormal <- xmlNode("LogNormal", attrs = c(id = paste0("LogNormalDistributionModel.", x$id),
                                              name = "distr"),
                       .children = list(xmlNode("parameter", 1,
                                                attrs = c(estimate = "false",
                                                          id = paste0("RealParameter.M", x$id),
                                                          name = "M")),
                                        xmlNode("parameter", 1.25,
                                                attrs = c(estimate = "false",
                                                          id = paste0("RealParameter.S", x$id),
                                                          lower = "0.0",
                                                          name = "S",
                                                          upper = "5.0"))
                                        ))
  
  xmlNode("distribution", taxonset, LogNormal,
          attrs = c(id = paste0(x$id, ".prior"),
                    spec = "beast.math.distributions.MRCAPrior",
                    tree = paste0("@Tree.t:", id)))
}