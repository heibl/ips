## This code is part of the ips package
## © C. Heibl 2014 (last update 2017-12-19)

#' @title Identification of MRCAs for Clades
#' @description Identify the most recent common ancestor (MRCA) nodes for one or
#'   more sets of taxa/tips.
#' @param phy An object of class \code{\link[ape]{phylo}}.
#' @param group A vector or list of vectors of mode \code{character} specifying
#'   the taxon set(s).
#' @param regex A logical, if \code{regex = TRUE}, taxon sets are matched to the
#'   tip labels as regular expressions of the form \code{"taxon1|taxon2"};
#'   otherwise strings will be matched exactly (see \code{\link{which}}).
#' @param stem Logical, ...
#' @param monophyletic Logical, ...
#' @return A vector of mode \code{"numeric"} containing node numbers.
#' @seealso \code{\link[ape]{mrca}}; \code{\link{descendants}} for the contrary operation to \code{noi}.
#' @examples 
#' # molecular phylogeny of Eurasian vipers:
#' # ---------------------------------------
#' #data(viperidae)	
#' #gen <- sapply(viperidae$tip.label, function(x) unlist(strsplit(x, "_"))[1])
#' #tax <- data.frame(genus = gen, species = viperidae$tip.label, row.names = NULL)
#' 
#' # group be a data frame
#' # ----------------------------------------------------
#' ## .. to be added ..
#' 
#' # group can be a list
#' # -------------------
#' #myclades <- split(tax$species, tax$genus)
#' #nds <- noi(viperidae, myclades)
#' #plot(viperidae)
#' #nodeInfo(nds)
#' 
#' # group might contain tip numbers
#' # -------------------------------
#' #group <- list(c(17, 22), c(13, 1))
#' #plot(viperidae)
#' #append2tips(phy, tips = unlist(group), pch = 19)
#' #nds <- noi(viperidae, myclades)
#' #nodeInfo(nds)
#' 
#' # the 'group' argument can also take regular expressions
#' # ------------------------------------------------------
#' #rex <- "aspis"
#' #node <- noi(viperidae, rex, regex = TRUE)
#' #plot.phylo(viperidae, tip.color = 0, edge.color = 0)
#' #box.clades(viperidae, nodes = node, col = "#D2A6A7", align = "all")
#' #plot.phylo.upon(viperidae)
#' #nodelabels(node = node, pch = 21, cex = 1.2, col = "red", bg = "#D2A6A7")
#' 
#' # if the 'group' argument is a list of elements of length 2,
#' # n = length(group) nodes of interest will be returned
#' # ----------------------------------------------------
#' #group <- list(
#' #    c("Vipera_berus", "Vipera_ursinii"),
#' #    c("Vipera_aspis_ssp._aspis", "Vipera_latastei"),
#' #    c("Vipera_ammodytes_ssp._ammodytes", 
#' #        "Vipera_ammodytes_ssp._montandoni"),
#' #    c("Macrovipera_lebetina", "Vipera_wagneri")
#' #)
#' #clades <- noi(viperidae, group)
#' #plot.phylo(viperidae, tip.color = 0, edge.color = 0)
#' #box.clades(viperidae, nodes = clades, col = c("#FFFFA5", "#D2A6A7",
#' #    "#A7D2A5", "#A5A6D2"), align = "all")
#' #plot.phylo.upon(viperidae)
#' @export

# 1. add option stemgroup
# 2. write example
# 3. check for nonvalid group members

noi <- function(phy, group, regex = FALSE, stem = FALSE, 
                monophyletic = FALSE){
  
  ## check phylogeny
  if ( !inherits(phy, "phylo") ) 
    stop("'phy' is not of class 'phylo'")  
	
	# test to filter out trees with 
  # non-consecutive tip labels
	# ------------------------------
  are.tips.consecutive(phy)
	
	# core function
	# -------------
	get.mrca <- function(phy, group, regex){
	  # apply regular expression matching
	  # ---------------------------------
    if ( regex ) {
      group <- phy$tip.label[grep(paste(group, collapse = "|"), phy$tip.label)]  
    }
    nn <- which(phy$tip.label %in% group)
    if ( length(nn) > 1) {
      repeat {
        nn <- sort(unique(phy$edge[phy$edge[, 2] %in% nn, 1]))
        gg <- phy$tip.label[descendants(phy, min(nn))]
        if ( all(group %in% gg) ) break
      } 
    }
    min(nn)
	} # end of get.mrca
	
	# turn 'group' to list
	if (!is.list(group)) group <- list(group)
  
  ## turn factor to character strings
  ## --------------------------------
  if (is.factor(group[[1]])){
    defactor <- function(x) levels(x)[x]
    group <- lapply(group, defactor)
  }
  
	## convert tip numbers to tip labels
	## ---------------------------------
	if (mode(group[[1]]) == "numeric" & !regex)
	  group <- lapply(group, function(x, phy) 
	    phy$tip.label[x], phy = phy)
  
 	# check tip labels
	# ----------------
	chk <- setdiff(unlist(group), phy$tip.label)
	if (length(chk) > 0 & !regex)						
	    stop(paste("tiplabels in 'group', but not in 'phy':\n", paste(chk, collapse = "\n")))
		
	nodes <- sapply(group, get.mrca, phy = phy, regex = regex)
  
  ## check for monophyly
  ## -------------------
  if (monophyletic){
    is.mp <- function(phy, node)
      phy$tip.label[descendants(phy, node)]
    tips <- lapply(nodes, is.mp, phy = phy)
    mp <- vector(length = length(tips))
    for ( i in seq_along(tips) ){
      mp[i] <- all(tips[[i]] %in% group[[i]])
    }
    nodes <- nodes[mp]
  }
  
  ## stem group nodes
  ## ----------------
  if ( stem ){
    nodes <- phy$edge[phy$edge[, 2] %in% nodes, 1]
  }
  
  return(nodes)
}
