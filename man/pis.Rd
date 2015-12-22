\name{pis}
\alias{pis}
\title{Number of Potentially-Informative Sites}

\description{
This function returns the number or positions of potentially-informative (parsimony-informative, phylogenetically-informative) sites in DNA sequence alignment.
}

\usage{
pis(x, what = "fraction", use.ambiguities = FALSE)
}

\arguments{
	
  \item{x}{ An object of class \code{\link{DNAbin}}.}
  
  \item{what}{Either of \code{"absolute"}, \code{"fraction"}, or \code{"index"}, which will return the absolute number, the relative number or the indeces of the potentially-informative sites.}
  
  \item{use.ambiguities}{\emph{Not yet available}.
  %Logical, currently only \code{FALSE} is implemented, which means that IUPAC ambiguity symbols will be treated as missing information (N).
  }
}

\value{
  Numeric (depending on \code{what}, the number, fraction, or indices of potentially-informative nucleotide sites).
}

\author{Christoph Heibl}

%\seealso{\code{\link{homogen.base.freq}}}

\examples{
# example data:
# -------------
data(ips.16S)	

# number of potentially-informative sites:
# ----------------------------------------
pis(ips.16S, what = "abs")
	
# proportion of potentially-informative sites:
# --------------------------------------------
pis(ips.16S, what = "frac")

# indeces of potentially-informative sites:
# --------------------------------------------
pis(ips.16S, what = "ind")
}
