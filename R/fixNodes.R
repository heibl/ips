## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-11-08)

#' @export

fixNodes <- function(phy){
  
  if (!inherits(phy, "phylo"))                             
    stop("object 'phy' is not of class 'phylo'")
  unfixed <- phy
  s <- 999999
  tip <- 1:Ntip(phy)
  internal <- Ntip(phy) + (1:Nnode(phy))
  pool <- union(tip, internal)
  root <- min(internal)
  is.internal <- function(x, maxtip){x > maxtip}
  phy$edge <- phy$edge + s
  
  ## fix root node
  ## -------------
  phy$edge[phy$edge[, 1] == min(internal) + s, 1] <- root
  pool <- setdiff(pool, root)
  internal <- setdiff(internal, root)
  
  ## fix non-root nodes
  ## ------------------
  for ( i in 1:nrow(phy$edge) ){
    n <- phy$edge[i, 2]
    if ( n <= Ntip(phy) + s ){
      phy$edge[i, 2] <- tip[1]
      tip <- tip[-1]
    } else {
      phy$edge[phy$edge[, 1] == n, 1] <- internal[1]
      phy$edge[i, 2] <- internal[1]
      internal <- internal[-1]
    }
  }
  
  ## fix tip- and nodelabels
  ## -----------------------
  trans <- data.frame(unfixed = unfixed$edge[, 2], 
                      fixed = phy$edge[, 2])
  ttrans <- trans[trans[, 2] <= Ntip(phy), ]
  ntrans <- trans[trans[, 2] > Ntip(phy), ]
  tid <- match(ttrans$unfixed, ttrans$fixed)
  phy$tip.label <- phy$tip.label[tid]
  nid <- match(ntrans$unfixed, ntrans$fixed)
  phy$node.label[-1] <- phy$node.label[-1][]
  
  # fix additional <node.label> elements
  # ------------------------------------
  sdtnames <- c("edge", "edge.length", "Nnode", "tip.label")
  if ( !all(names(phy) %in% sdtnames) ){
    nls <- which(!names(phy) %in% sdtnames)
    for ( i in nls )
      #phy[[i]] <- phy[[i]][id]
      phy[[i]][-1] <- phy[[i]][-1][nid]
  }  
  phy
}