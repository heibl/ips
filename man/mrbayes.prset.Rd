\name{mrbayes.prset}
\alias{mrbayes.prset}
\title{Set Priors for MrBayes}
\description{
Set prior distributions for \code{\link{mrbayes}}.
}
\usage{
mrbayes.prset(...)
}
\arguments{
  \item{\dots}{
arguments in \code{tag = value} form, or a list of tagged values. The tags must come from the names of prior distribution parameters described in the ‘Prior Distribution Parameters’ section.
}
}
%\details{}
\value{
a list of length zero (see 'Note')
%a list containing a subset (including the empty and the full set) of model parameters.
}
\section{Prior Distribution Parameters}{
\describe{
 \item{traitiopr}{}
 \item{revmatpr}{}
 \item{aamodelpr}{}
 \item{aarevmatpr}{}
 \item{omegapr}{}
 \item{ny98omega1pr}{}
 \item{ny98omega3pr}{}
 \item{m3omegapr}{}
 \item{codoncatfreqs}{}
 \item{statefreqpr}{}
 \item{shapepr}{}
 \item{ratecorrpr}{}
 \item{pinvarpr}{}
 \item{covswitchpr}{}
 \item{symdirihyperpr}{}
 \item{topologypr}{}
 \item{brlenspr}{}
 \item{clockvarpr}{}
 \item{igrvarpr}{}
  }
  }
\references{
J.P. Huelsenbeck & Ronquist F. 2001. MrBayes: Bayesian inference of phylogenetic trees. \emph{Bioinformatics} \bold{17}: 754-755.

Ronquist F. & J.P. Huelsenbeck. 2003. MrBayes 3: Bayesian phylogenetic inference under mixed models. \emph{Biometrics} \bold{19}: 1572-1574.
	
MrBayes website: \url{https://mrbayes.sourceforge.net/}.
}
\author{Christoph Heibl}

\note{This function currently returns an empty set of prior distribution parameters, i.e., you cannot change the MrBayes default parameters. }
\seealso{
\code{\link{mrbayes.lset}} to set model parameters, \code{\link{mrbayes.mcmc}} to set parameters of the Markov chain Monte Carlo (MCMC), and \code{\link{mrbayes}} to run MrBayes locally or prepare input files for a computer cluster.
}
\examples{
mrbayes.prset()
}