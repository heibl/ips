\name{tipHeights}
\alias{tipHeights}
\title{
Tip Heights in a Phylogenetic Tree
}
\description{
For each tip (leave, terminal node) in the phylogenetic tree the edge lengths (branch lengths) from root to tip, be it units of time or divergence, is summed up.
}
\usage{
tipHeights(phy)
}
\arguments{
  \item{phy}{an object of class \code{\link[ape:read.tree]{phylo}}.}
}

\value{
a numeric vector with distances from root to tip for each tip in the phylogenetic tree.
}
\author{
Christoph Heibl
}
\seealso{
\code{\link[ape]{branching.times}}
}
%\examples{}

