\name{prank}
\alias{prank}

\title{PRANK}

\description{
DNA sequence Alignment Using the program PRANK.}

\usage{
prank(x, outfile, guidetree = NULL, gaprate = 0.025, 
    gapext = 0.75, path)
}

\arguments{
	
  \item{x}{an object of class \code{DNAbin}.}
  
  \item{outfile}{a character string giving a name for the output file.}
  
  \item{guidetree}{an object of class \code{phylo} to be used as guidetree in alignment.}
  
  \item{gaprate}{numeric giving the gap opening rate; defaults to 0.025.}
  
  \item{gapext}{numeric giving the gap extension penalty; defaults to 0.75.}
  
  \item{path}{a character string indicating the path to the PRANK executable.}
  
}

\value{
  \code{matrix} of class \code{"DNAbin"}
}

\references{
	\url{http://www.ebi.ac.uk/goldman-srv/prank/prank/}
}

\note{\code{prank} was last updated and tested to work with PRANK v. 120814 on Windows XP. If you have problems getting the function to work with a newer version of PRANK, contact the package maintainer.}

\author{Christoph Heibl}

\seealso{
	\code{\link{read.fas}} to import DNA sequences;
	\code{\link{mafft}} for another alignment algorithm;
	\code{\link{gblocks}} and \code{\link{aliscore}} for alignment cleaning. 
}
