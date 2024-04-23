\name{read.beast}
\alias{read.mrbayes}
\alias{read.beast}
\alias{read.starbeast}
\title{Read Bayesian Trees}

\description{
Thess functions parse chronograms in NEXUS format as produced by TreeAnnotator or output by MrBayes.}

\usage{
read.mrbayes(file, digits = NULL)

read.beast(file, digits = NULL)

read.starbeast(file)
}

\arguments{
	
  \item{file}{A character string giving the input file, which must be a TreeAnnotator-generated chronogram in NEXUS format.}
  
  \item{digits}{NULL or integer, if \code{!is.null(digits)} values are rounded to the given integer.}
}

\value{
  An object of class \code{phylo}.
}

\note{
\code{read.starbeast} currently parses only skalars and ranges; node statistics with more than two values will be deleted and a warning message will be issued. Future version of \code{read.starbeast} will hopefully be able to append list or data frames to \code{phylo} objects. If you have any opinion or wishes regarding the question of how this exactly should be managed, send me a message.

\code{read.mrbayes} is intended to read single trees with annotated nodes. For reading tree samples from the posterior distribution there is a collection of perl script written by Johan Nylander (\url{https://github.com/nylander/Burntrees}).
}

\references{
TreeAnnotator: \url{http://beast.community/treeannotator}

Metacomments in NEXUS: \url{https://code.google.com/archive/p/beast-mcmc/wikis/NexusMetacommentFormat.wiki}	
}

\author{Christoph Heibl}

\seealso{\code{\link{read.beast.table}} to extract internal node data from NEXUS file, \code{\link{rbeauti}} to create XML input for BEAST. \code{HPDbars} for plotting highest posterior densities on phylogenies has been moved to package \bold{viper}.
}
