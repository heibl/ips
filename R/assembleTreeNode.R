## This code is part of the ips package
## © C. Heibl 2014 (last update 2020-03-10)

#' @export

assembleTreeNode <- function(i, x){
  
  ## If trees are linked, this node is only required for the first partition
  ## ----------------------------------------------------------------------
  if (i > 1 & x$link.trees) return(NULL)
  
  ## Prepare id attributes: they equal taxon_set_id, except the tree id for
  ## linked trees
  ## ------------
  general_id <- x$partition[i]
  # tree_id <- ifelse(x$link.trees, "trees", taxon_set_id)
  tree_id <- ifelse(x$link.trees, "trees", general_id)
  
  if (x$tree %in% c("Yule", "Calibrated Yule", "Birth Death")){
    
    Tree.t <- list(
      xmlNode("taxonset", 
              attrs = c(id = paste0("TaxonSet.", general_id), spec = "TaxonSet"),
              .children = list(xmlNode("alignment",
                                       attrs = c(idref = general_id))))
    )
  
  # } else if (x$tree == "xxxx"){
  #   
  #   Tree.t <- list(
  #     xmlNode("trait", 
  #             attrs = c(id = paste0("dateTrait.t:", general_id), 
  #                       spec = "beast.evolution.tree.TraitSet",
  #                       traitname = "date-backward",
  #                       value = x$tip.dates),
  #             .children = list(xmlNode("taxa", 
  #                                      attrs = c(id = paste0("TaxonSet.", general_id),
  #                                                spec = "TaxonSet"),
  #                                      .children = list(xmlNode("alignment",
  #                                                               attrs = c(idref = general_id)))))),
  #     xmlNode("taxonset", attrs = c(idref = paste0("TaxonSet.", general_id)))
  #   )
    
  } else if (x$tree == "Fossilized Birth Death"){
    
    ## Fossilized Birth Death prior
    ## ----------------------------
    Tree.t <- list(
      xmlNode("trait", 
              attrs = c(id = paste0("dateTrait.t:", general_id), 
                        spec = "beast.evolution.tree.TraitSet",
                        traitname = "date-backward",
                        value = x$tip.dates),
              .children = list(xmlNode("taxa", 
                                       attrs = c(id = paste0("TaxonSet.", general_id),
                                                         spec = "TaxonSet"),
                                       .children = list(xmlNode("alignment",
                                                                attrs = c(idref = general_id)))))),
      xmlNode("taxonset", attrs = c(idref = paste0("TaxonSet.", general_id)))
    )
    
  } else {
    stop("incorrect specification of tree prior")
  }
  
  Tree.t <- xmlNode("tree", 
                    attrs = c(id = paste0("Tree.t:", general_id),
                              spec = "beast.evolution.tree.Tree",
                              name = "stateNode"),
                    .children = Tree.t)
  
  Tree.t
}
