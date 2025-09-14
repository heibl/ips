\name{rbeauti}
\alias{rbeauti}

\title{XML Input Files for BEAST}

\description{
This function is intended to prepare XML files for BEAST with R. BEAST uses an MCMC approach to estimate rooted phylogenies from molecular data (Drummond & Rambaut, 2007). %See \code{note} section for currently supported models.
}

\usage{
rbeauti(..., file, template = "standard", taxonset)
}

\arguments{

  \item{...}{one or more object(s) of class \code{\link[ape]{DNAbin}}.}
	
  \item{file}{A connection, or a character string naming the file to write to. If left empty the XML tree will be printed to the screen (see Examples).}
  \item{template}{\emph{Currently unused.}}
 
  \item{taxonset}{A list containing one or more taxon sets.}
%  
%  \item{monophyly}{A vector indicating monophyly constraints for the taxon sets declared with \code{taxonset}.}
%  
%  \item{tmrcaCons}{A list, containing the prior distribution(s) for age constraints of internal nodes (which must be grouped by \code{taxonset}).}
%  
%  \item{startingTree}{Either "random" or "upgma", or an object of class \code{"phylo"} to be used as a starting tree.}
%  
%  \item{specModel}{A character string indicating an evolutionary model to construct a prior distribution of node heights (tree prior). Currently implemented are the \bold{Yule Model} (\code{"yule"}), the \bold{Birth-Death Model} (\code{"birthDeath"}), and the \bold{Coalescent model with constant size} (choosen with any other string).}
%  
%  \item{clock}{A character string, either \code{"strict"} or \code{"lognormal"} to choose between the strict clock and the uncorrelated lognormal relaxed clock model. The exponential rates relaxed clock is currently not supported.}
%  
%  \item{ngen}{A character string, the number of generations to run the MCMC.}
%  
%  \item{samplefreq}{A character string, the intervals between sampling the MCMC.}
%  
%  \item{logSubs}{A logical, indicating if trees with branch lengths expressed as substitution should be logged.}
%  
%  \item{nodata}{A logical, indicating if BEAST should be run without data (see details).} 
}
%
\details{
\code{rbeauti} has been completely rewritten to work with \bold{BEAST 2}. Currently \code{rbeauti} offers few options, because the idea is not to create ready-to-use XML file. That can be done conveniently with \bold{BEAUti} (the BEAST package's genuine XML generator). Instead, \code{rbeauti} is intended to make the definition of large numbers of taxon sets easy. The creation of taxon sets can be done via R scripts and the resulting XML files can be further modified with BEAUti.
}
%
\references{
The BEAST 2 website: \url{https://www.beast2.org/}	

Drummond, A.J. & A. Rambaut. 2007. BEAST: Bayesian evolutionary analysis by sampling trees. \emph{BMC Evolutionary Biology} \bold{7}: 240.
}
%
\author{Christoph Heibl}

%\note{}

\seealso{
  \code{\link{read.beast}}, \code{\link{read.beast.table}}
}
\examples{
data(ips.16S)

## define taxon sets
spec <- rownames(ips.16S)
ingroup <- spec[grep("Ips|Orthomotomicus", spec)]
outgroup <- spec[grep("Pityogenes", spec)]
ts <- list(list(id = "ingroup", taxon = ingroup),
           list(id = "outgroup", taxon = outgroup))

## print XML file to screen
rbeauti(ips.16S, taxonset = ts)
}
