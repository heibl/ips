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
  if (any(info)){
    
    ## one column perfectly matching
    tip_col <- names(info)[info]
    stopifnot(length(tip_col) == 1)
    message("tip labels are taken from column '", tip_col, "'")
    
  } else {
    
    ## no perfect match between tip labels and any column of 'data'
    info <- apply(data, 2, function(a, b) which(a %in% b), b = phy$tip.label)
    n <- sapply(info, length)
    tip_col <- names(n)[which.max(n)]
    not_in_data <- setdiff(phy$tip.label, data[[tip_col]])
    n_miss <- length(not_in_data)
    not_in_phy <- setdiff(data[[tip_col]], phy$tip.label)
    info <- data.frame(phy = not_in_data, data = not_in_phy)
    id <- sapply(info$phy, agrep, x = info$data)
    info$data <- info$data[id]
    info <- paste0(" - '", info$phy, "' (phy) vs. '", info$data, "' (data)")
    info <- paste(info, collapse = "\n")
    message("tip labels seem to be stored in column '", tip_col, "', but ", 
            ifelse(n_miss == 1, "1 does", paste(n_miss, "do")), " not match:")
    message(info)
    stop("'data' is not congruent with tiplabels of 'phy'")
  }
  
  ## Calculate number of accessions per species
  ## ------------------------------------------
  comb_or_not <- table(data[, tip_col])
  
  data <- split(data[, names(data) != tip_col], data[, tip_col])
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