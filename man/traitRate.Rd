\name{traitRate}
\alias{traitRate}
\title{Trait-Dependent Shifts in Molecular Rate}
\description{
Detection of trait-dependent shifts in the rate of molecular evolution with \bold{traitRate} (Mayrose & Otto, 2011).
}
\usage{
traitRate(phy, seq, x, mainType = "Optimize_Model", 
   n, charModelParam1 = 0.5, charModelParam2 = 1, 
   gammaParam = 0.5, seqModelParam1 = 2,
   exec = "/Applications/traitRate-1.1/programs/traitRate")
}
\arguments{
  \item{phy}{a ultrametric phylogenetic tree of class \code{\link{phylo}.}
}
  \item{seq}{a multiple sequence alignment of class \code{\link{DNAbin}}.
}
  \item{x}{data frame containing a binary character in the first column.
}
  \item{mainType}{character string giving the type of analysis; two choices are possible: \code{"Optimize_Model"} will produce MLE of parameters and \code{"runTraitBootstrap"} will perform a parametric bootstrap analysis.}
  \item{n}{numeric, the number of bootstrap replicates. Will be ignored if \code{mainType = "Optimize_Model"}.}
  \item{charModelParam1}{numeric, giving an initial value for the rate of transitions of character state 0 to 1.}
  \item{charModelParam2}{numeric, giving an initial value for the rate of transitions of character state 1 to 0.}
  \item{gammaParam}{numeric, giving an initial value for the \bold{alpha} parameter of the model of sequence evolution.}
  \item{seqModelParam1}{numeric, giving an initial value for the \bold{kappa} parameter of the model of sequence evolution.}
  \item{exec}{character string giving the path to the program directory.
}
}
%\details{}
\value{
Currently none, but look for the output files in the \code{'RESULTS'} subdirectory in the current working directory.
}
\references{
Mayrose, I. & S.P. Otto. 2011. A likelihood method for detecting trait-dependent shifts in the rate of molecular evolution. \emph{Mol. Biol. Evol.} \bold{28}: 759-770
}
\author{
Christoph Heibl
}
\note{
This function is under development!
}

\seealso{
\code{\link{read.tree}} for reading phylogenetic trees, \code{\link{read.fas}} for reading multiple sequence alignments in FASTA format.
}
%\examples{}
