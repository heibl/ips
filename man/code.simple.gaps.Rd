\name{code.simple.gaps}
\alias{code.simple.gaps}

\title{Simple Gap/Indel Coding}

\description{
\code{code.simple.gaps} takes an aligned DNA sequence matrix and codes the simple gaps, i.e. gaps that do not overlap with other gaps. The gapped positions are excluded from the matrix and the coded gap characters are appended to the matrix.
}

\usage{
code.simple.gaps(x, append = TRUE)
}

\arguments{
  \item{x}{An object of class \code{\link[ape]{DNAbin}}.}
  
  \item{append}{Logical.}
}

\value{
An object of class \code{\link[ape]{DNAbin}.}
}
 
\references{
Simmons, M.P. & H. Ochoterena. 2000. Gaps as characters in sequence-based phylogenetic analyses. \emph{Systematic Biology} \bold{49(2)}: 369--381.	
}

\author{Christoph Heibl}

\seealso{
	\code{\link{deleteGaps}}, \code{\link{deleteEmptyCells}}, \code{\link{trimEnds}}
	}
