\name{deleteGaps}
\Rdversion{1.1}
\alias{deleteGaps}
\alias{delete.gaps}

\title{Remove Gap/Indel Positions from Alignment}

\description{
This function removes indel positions (or gaps) from a DNA sequence alignment. For faster execution, \code{deleteGaps} handles sequences in \bold{ape}'s bit-level coding scheme.
}

\usage{
deleteGaps(x, nmax = nrow(x) - 4)
}

\arguments{
  \item{x}{an object of class \code{\link[ape]{DNAbin}}.}
  \item{nmax}{an integer number between \code{0} and \code{nrow(x) - 2}, which gives the maximum number of
  gap characters ("-") that will be tolerated at any given alignment position (column).}
}
\details{
The default, \code{nmax = nrow(x) - 4}, removes all those positions from the alignment, which contain at least four non-gap characters, which is the minimum number of sequences needed to produce a non-trivial unrooted topology. All gaps will be excluded by selecting \code{nmax = 0} and half of all gaps with \code{nmax = nrow(x) / 2}.  

In contrast, \code{\link[ape]{del.gaps}} removes all gap characters from the alignment, so most probably the result will not be a set of sequences of equal length and the matrix will be coerced to a list.
}
\value{
an object of class \code{\link[ape]{DNAbin}.}
}

\author{Christoph Heibl}

\seealso{
\code{\link{code.simple.gaps}} for coding of simple gaps, \code{\link[ape]{del.gaps}} for removal of all gap symbols from an alignment, \code{\link{gblocks}} and \code{\link{aliscore}} for more sophisticated  methods of cleaning/masking alignments.
}