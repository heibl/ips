## This code is part of the ips package
## © C. Heibl 2014 (last update 2017-12-20)

#' @title Identification of Stem-Lineage-Edges and MRCAs
#' @description \code{noi} (\strong{n}ode \strong{o}f \strong{i}nterest)
#'   identifies the most recent common ancestor (MRCA) and \code{eoi}
#'   (\strong{e}dge \strong{o}f \strong{i}nterest) its subtending stem-lineage
#'   edge of one or more sets of taxa/tips.
#' @param phy An object of class \code{\link[ape]{phylo}}.
#' @param node A vector of mode \code{"numeric"} giving the nodes numbers of the nodes whose subtending stem-lineages will be identified.
#' @param group A vector or list of vectors of mode \code{character} specifying
#'   the taxon set(s). Will be ignored if \code{node} is given.
#' @param regex A logical, if \code{regex = TRUE}, taxon sets are matched to the
#'   tip labels as regular expressions of the form \code{"taxon1|taxon2"};
#'   otherwise strings will be matched exactly (see \code{\link{which}}).
#' @param stem Logical, ...
#' @param monophyletic Logical, ...
#' @return A vector of mode \code{"numeric"} containing node numbers.
#' @seealso \code{\link[ape]{mrca}}; \code{\link{descendants}} for the contrary
#'   operation to \code{noi}.
#' @examples
#' # molecular phylogeny of Ips bark beetles
#' # ---------------------------------------
#' data(ips.tree)
#' ips.tree <- ladderize(ips.tree)
#' ips.tree <- fixNodes(ips.tree)
#' clade1 <- descendants(ips.tree, 44, labels = TRUE)
#' mrca <- noi(ips.tree, clade1)
#' stem_lineage <- eoi(ips.tree, mrca)
#' ecol <- rep("black", Nedge(ips.tree))
#' ecol[stem_lineage] <- "red"
#' plot(ips.tree, no.margin = TRUE, edge.color = ecol)
#' nodelabels(node = mrca, pch = 22, col = "blue")
#' #gen <- sapply(viperidae$tip.label, function(x) unlist(strsplit(x, "_"))[1])
#' #tax <- data.frame(genus = gen, species = viperidae$tip.label, row.names = NULL)
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
#' @name oi

NULL