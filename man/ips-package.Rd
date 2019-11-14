\name{ips-package}
\alias{ips-package}
\alias{ips}
\docType{package}
\title{Interfaces to Phylogenetic Software}
\description{
This package presents a set of functions that were formerly included in the \emph{phyloch} package and which wrap popular phylogenetic software for sequence alignment, masking of sequence alignments, and estimation of phylogenies and ancestral character states.
}
\details{
\tabular{ll}{
 Package: \tab ips\cr
 Type: \tab Package\cr
Version: \tab 0.0.12\cr
Date: \tab 2019-11-14\cr
 License: \tab GPL (>= 2)\cr
}
There are several functions for reading and writing DNA sequences in FASTA, PHYLIP, and NEXUS format: \code{\link{read.fas}}, \code{\link{read.phy}}, \code{\link{read.nex}}, \code{\link{write.fas}}, \code{\link{write.phy}}, and \code{\link{write.nex}}.
Some functions are available for integrating BEAST with R. XML input files for BEAST can be generated with \code{\link{rbeauti}}. Two functions are designed to read TreeAnnotator output: \code{\link{read.beast}} will render an object of class \code{phylo} with additional node statistics appended as list elements. These additional node statistics will be lost be the subsequent use of \code{\link{ladderize}} or \code{\link{rotate}} (or similar functions that change the ordering of internal nodes).\code{\link{read.beast.table}} also parses the TreeAnnotator output, but returns a matrix of node statistics. 
This package itself does not implement techniques for phylogenetic analyses, but provides a series of wrappers for commonly used software packages. Sequence alignment can be done with the \code{\link{mafft}} and \code{\link{prank}}; cleaning of sequences with \code{\link{gblocks}} and \code{\link{aliscore}}. The function \code{\link{raxml}} and \code{\link{mrbayes}} are intended for phylogenetic tree search. Running \code{\link{mrbayes}} with argument \code{run = FALSE} can be used to create MrBayes-executable NEXUS files. Finally,  wrappers is provided for \code{Multistate} in the \code{BayesTraits} package (see \code{\link{multistateML}} and \code{\link{multistateMCMC}}).
Several plotting functions (\code{HPDbars}, \code{clade.bars}, \code{box.clades}, \code{box.tips}, \code{tip.color}, \code{edge.color} have been moved to the \bold{viper} package.
}
\author{Natalie Cusimano, Christoph Heibl, Franz-Sebastian Krah,
Maintainer: Christoph Heibl (\email{christoph.heibl@gmx.net})
}
\keyword{package}
\seealso{
\code{\link[ape:ape-package]{ape}}
}
