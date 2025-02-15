## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2025-02-15)

#' @title Standard Node Numbering in Phylo Objects
#' @description (Re-)establishes the standard numbering of terminal and internal
#'   nodes in phylogenies represented as objects of class \code{\link{phylo}}.
#' @param phy An object of class \code{\link{phylo}}.
#' @details 
#' 
#' Internal and terminal nodes of \code{\link{phylo}} objects have a standard,
#' "canonical" node numbering. It consists of the following rules:
#' \enumerate{
#'  \item Number tip nodes from left/bottom to right/top by 1 to \code{ntip()}
#'  \item Number internal nodes by preorder traversal (see Examples) consecutively
#'  from \code{ntip() + 1} to \code{nnode() + 1}
#' }
#'  
#' When reading phylogenetic trees from a NEXUS file that contains a
#' \code{translate} section, it can happen that the terminal nodes (tips,
#' leaves) of the corresponding \code{phylo} object are not numbered
#' consecutively, which can be a problem in some downstream applications. You
#' can use \code{fixNodes} to get the correct order of terminal node numbers.
#' 
#' \code{fixNodes} is also intended to re-establish the standard numbering of
#' internal nodes and reorder all node value elements (e.g. node.label,
#' posterior, ...) if a \code{\link[ape]{phylo}} object has been modified by
#' either \code{\link[ape]{root}}, \code{\link[ape]{ladderize}}, or
#' \code{\link[ape]{rotate}}.
#' @return An object of class \code{\link{phylo}}.
#' @seealso \code{\link[ape]{read.tree}}, \code{\link{read.nexus}},
#'   \code{\link{read.beast}} for reading trees in NEWICK and NEXUS format;
#'   \code{\link[ape]{ladderize}} and \code{\link[ape]{rotate}} for tree
#'   manipulation; \code{node.support} for plotting node support values has been
#'   moved to package \bold{viper}. %\code{\link{node.trans}} for handling
#'   non-standard node elements in class \code{\link[ape]{phylo}}.
#' @examples
#' set.seed(9999)
#' 
#' # Create random topology with random node numbering
#' phy1 <- rtopology(10)
#' 
#' # Impose 'canonical' node numbering
#' phy2 <- fixNodes(phy1)
#' 
#' # Compare node numbering schemes
#' op = par(no.readonly = TRUE)
#' par(mfcol = c(1, 2), mar = c(0, 0, 3, 0))
#' plot(phy1, use.edge.length = FALSE, show.tip.label = FALSE)
#' nodelabels(); tiplabels()
#' title("random node numbers")
#' plot(phy2, use.edge.length = FALSE, show.tip.label = FALSE)
#' nodelabels(); tiplabels()
#' title("preorder traversal of internal nodes")
#' par(op)
#' 
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