\name{descendants}
\alias{descendants}

\title{Descendants of an Internal Node in a Phylogeny}

\description{
For any given internal node of a phylogeny, the function returns a vector containing the node numbers descending from that node.
}

\usage{
descendants(phy, node, type = "t", ignore.tip = TRUE, labels)
}

\arguments{
  \item{phy}{an object of class \code{\link[ape]{phylo}}.}
  
  \item{node}{an integer giving the number of the internal node.}
  
  \item{type}{a character string, may be \code{"terminal"}, \code{"internal"}, \code{"both"}, or any unambiguous abbreviation of these.}
  
  \item{ignore.tip}{logical, if \code{ignore.tip = FALSE}, the function will issue an error when \code{node} is not internal, otherwise the number of the corresponding terminal node will be returned.}
  \item{labels}{logical, determines if node labels are returned instead of node number, currently ignored unless \code{type = "t"}.}
}

\value{
  A vector containing terminal node numbers or tip labels.
}

\author{Christoph Heibl}

\seealso{\code{\link{sister}}, \code{\link{noi}}}

\examples{
 # generate a random tree with 12 terminal and 11 internal nodes:
 tree <- rtree(12)
 
 # get the descendants of internal node 15:
 x <- descendants(tree, 15)
}