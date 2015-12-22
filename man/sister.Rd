\name{sister}
\alias{sister}

\title{Identification of Sister Nodes and Clades}

\description{For any given internal node in a phylogeny, this function returns the sister clade}

\usage{
sister(phy, node, type = "terminal")
}

\arguments{
  \item{phy}{an object of class \code{phylo}.}
  
  \item{node}{an integer giving the number of the node.}
  
  \item{type}{a character string, may be \code{"terminal"}, \code{"internal"}, \code{"both"}, or any unambiguous abbreviation of these.}
}

\value{
  A vector containing node numbers.
}

\author{Christoph Heibl}

\seealso{\code{\link{descendants}}, \code{\link{noi}}}

\examples{
# data
# ----
#data(vipera)
#node <- noi(vipera.ml, c("Vipera_latastei", "Vipera_aspis"))

# get the sister node of 'node'
# -----------------------------
#(sn <- sister(vipera.ml, node, type = "internal"))

# ... and the sister nodes's tips (i.e. the sister clade)
# -------------------------------------------------------
#(sc <- sister(vipera.ml, node, type = "terminal"))

# results graphically:
# --------------------
#plot(vipera.ml)
#nodelabels(text = "15", node = 15, frame = NULL, col = "white", bg = "black")
#nodelabels(text = sn, node = sn, col = "black", bg = "red")
#append2tiplabel(tree, tips = sc, pch = 21, offset = strwidth("."))
}