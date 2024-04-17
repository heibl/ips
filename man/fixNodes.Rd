\name{fixNodes}
\Rdversion{1.1}
\alias{fixNodes}

\title{
Standard Node Numbering in Phylo Objects
}
\description{
The function (re-)establishes the standard numbering of terminal and internal nodes in phylogenies represented as objects of class \code{\link{phylo}}.
}

\usage{
fixNodes(phy)
}
 
\arguments{
  \item{phy}{An object of class \code{\link{phylo}}. 
	}
}

\details{
When reading phylogenetic trees from a NEXUS file that contains a \code{translate} section, it can happen that the terminal nodes (tips, leaves) of the corresponding \code{phylo} object are not numbered consecutively, which can be a problem in some downstream applications. You can use \code{fixNodes} to get the correct order of terminal node numbers.

\code{fixNodes} is also intended to re-establish the standard numbering of internal nodes and reorder all node value elements (e.g. node.label, posterior, ...) if a \code{\link[ape]{phylo}} object has been modified by either \code{\link[ape]{root}}, \code{\link[ape]{ladderize}}, or \code{\link[ape]{rotate}}.
}

\value{An object of class \code{\link{phylo}}.}

\author{Christoph Heibl}

\note{
  \code{fixNodes} has been completely rewritten for \bold{ips} version 1.0-0. It should now run absolutely stable and is much quicker. Nevertheless, I recommend checking carefully the results of \code{fixNodes}, until the function has been tested by a number of users. Then this comment will be removed.
}

\seealso{\code{\link[ape]{read.tree}}, 
	\code{\link{read.nexus}}, \code{\link{read.beast}} 
  for reading trees in NEWICK and NEXUS format; \code{\link[ape]{ladderize}} and \code{\link[ape]{rotate}} for tree manipulation; \code{node.support} for plotting node support values has been moved to package \bold{viper}. %\code{\link{node.trans}} for handling non-standard node elements in class \code{\link[ape]{phylo}}.
}

%\examples{
%set.seed(9999)
%phy <- rtree(20)
%is.tiporder.consecutive(phy)
%phy <- ladderize(phy)
%is.tiporder.consecutive(phy)
%phy <- fixNodes(phy)
%is.tiporder.consecutive(phy)
%}
