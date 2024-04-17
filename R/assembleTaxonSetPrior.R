

assembleTaxonSetPrior <- function(taxonset, setname, id, i, link.trees){
  
  ## Attributes of <taxon> change if more the same MRCA prior 
  ## when there is more than one partition
  
  if (link.trees){
    
    if (i == 1){
      
      taxonSetPrior <- sapply(taxonset, 
                              function(x) xmlNode("taxon", attrs = c(id = x, spec = "Taxon")))
      taxonSetPrior <- xmlNode("taxonset", attrs = c(id = paste(setname, "trees", sep = "."), 
                                                     spec = "TaxonSet"),
                               .children = taxonSetPrior)
      
      taxonSetPrior <- xmlNode("distribution", 
                               attrs = c(id = paste(setname,"trees.prior", sep = "."),
                                         spec = "beast.math.distributions.MRCAPrior",
                                         monophyletic="true",
                                         tree = "@Tree.t:trees"),
                               .children = list(taxonSetPrior))
    } else {
      taxonSetPrior <- NULL
    }
    
  } else {
    
    ## Unlinked trees
    ## --------------
    if (i == 1){
      taxonSetPrior <- sapply(taxonset, 
                              function(x) xmlNode("taxon", attrs = c(id = x, spec = "Taxon")))
    } else {
      taxonSetPrior <- sapply(taxonset, 
                              function(x) xmlNode("taxon", attrs = c(idref = x)))
    }
    taxonSetPrior <- xmlNode("taxonset", attrs = c(id = paste(setname, id[i], sep = "."), 
                                                   spec = "TaxonSet"),
                             .children = taxonSetPrior)
    
    taxonSetPrior <- xmlNode("distribution", 
                             attrs = c(id = paste(setname, id[i], "prior", sep = "."),
                                       spec = "beast.math.distributions.MRCAPrior",
                                       monophyletic="true",
                                       tree = paste0("@Tree.t:", id[i])),
                             .children = list(taxonSetPrior))
    
  }
  taxonSetPrior
}