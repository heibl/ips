## This code is part of the ips package
## © C. Heibl 2014 (last update 2017-12-20)

#' @rdname oi
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
