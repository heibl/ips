## This code is part of the ips package
## © C. Heibl 2018 (last update 2020-12-16)

#' @title Graft Polytomies on Tips of Phylogeny
#' @description Graft polytomies on the tips of a class \code{phylo} object.
#' @param phy An object of class \code{\link{phylo}}.
#' @param data A data frame containing two columns. The entries of one column
#'   must be identical to the tip labels of the phylogeny; the other column
#'   contains the new tip labels. The column are matched
#'   automatically.
#' @param brlen A numeric giving the branch lengths for the polytomies.
#' @param annotate Logical, if \code{TRUE}, the former tip labels will be turned
#'   into node labels. Note, that this will overwrite existing node labels.
#' @return An object of class \code{\link{phylo}} with \code{nrow(data)} tips.
#' @seealso \code{\link{forceEqualTipHeights}}
#' @examples
#' data(ips.tree)
#' ## Simulate varying number of intraspecific observations
#' s <- sapply(1:Ntip(ips.tree), sample, x = 1:3, size = 1)
#' x <- rep(ips.tree$tip.label, times = s)
#' x <- data.frame(x, paste0(x, unlist(lapply(s, function(z) 1:z))))
#' ## Create polytomies
#' tre <- combMyTree(ips.tree, x)
#' plot(tre, no.margin = TRUE, cex =.5)
#' @importFrom ape bind.tree compute.brlen read.tree
#' @export

combMyTree <- function(phy, data, brlen = 0, annotate = FALSE){
  
  ## Do some tests
  ## -------------
  if (!inherits(phy, "phylo")) stop("'phy' is not of class 'phylo'")
  
  info <- apply(data, 2, setequal, y = phy$tip.label)
  if (!any(info)) stop("'data' is not congruent with tiplabels of 'phy'")
  
  ## Calculate number of accessions per species
  ## ------------------------------------------
  comb_or_not <- table(data[, info])
  
  ## Case 1: Only 1 accession, simply replace
  ##         species name by accession name
  ## --------------------------------------
  # if (any(comb_or_not == 1)){
  #   one_sample <- names(comb_or_not)[comb_or_not == 1]
  #   one_sample <- phy$tip.label[phy$tip.label %in% one_sample]
  #   phy$tip.label[phy$tip.label %in% one_sample] <- data[match(one_sample, data[, info]), !info]
  # }
  
  ## Case 2: More than 1 accession, put combs on tips
  ## ------------------------------------------------
  # data <- data[data[, info]  %in% names(comb_or_not)[comb_or_not > 1], ]
  data <- split(data[, !info], data[, info])
  makeComb <- function(z){
    z <- read.tree(text = paste0("(", paste(z, collapse = ","), ");"))
    compute.brlen(z, brlen) ## set branch lengths to br.len
  }
  combs <- lapply(data, makeComb)
  for (i in seq_along(combs)){
    phy <- bind.tree(phy, combs[[i]], which(phy$tip.label == names(combs)[i]))
  }
  phy <- fixNodes(phy)
  
  ## Annotate previous tips
  ## ----------------------
  for (i in seq_along(data)){
    phy$node.label[noi(phy, data[[i]]) - Ntip(phy)] <- names(data)[i]
  }
  
  phy
}