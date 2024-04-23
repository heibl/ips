\name{mrbayes.lset}
\alias{mrbayes.lset}
\title{Model Settings for MrBayes}
\description{
Set model parameters for \code{\link{mrbayes}}.
}
\usage{
mrbayes.lset(..., partition)
}
\arguments{
  \item{\dots}{
arguments in \code{tag = value} form, or a list of tagged values. The tags must come from the names of model parameters described in the ‘Model Parameters’ section.
}
  \item{partition}{
a character string giving the labelling for a partion.
}
}
%\details{}
\value{
a list containing a subset (including the empty and the full set) of model parameters.
}
\section{Model Parameters}{
\describe{
  \item{nucmodel}{
\code{"4by4"}, \code{"doublet"}, \code{"codon"}, or \code{"protein"}.
  }
  \item{nst}{
\code{1}, \code{2}, \code{6}, or \code{"mixed"}.
  }
  \item{code}{
\code{"universal"}, \code{"vertmt"}, \code{"mycoplasma"}, \code{"yeast"}, \code{"ciliates"}, or \code{"metmt"}. 
  }
  \item{ploidy}{
\code{"haploid"}, \code{"diploid"},  or \code{"zlinked"}.
  }
    \item{rates}{
\code{"equal"}, \code{"gamma"}, \code{"propinv"},\code{"invgamma"}, or \code{"adgamma"}.
  }
  \item{ngammacat}{
\code{1-24} 
  }
  \item{nbetacat}{
\code{1-24}
  }
  \item{omegavar}{
\code{"equal"},  \code{"ny98"},  or \code{"m3"}.
  }
  \item{covarion}{
\code{"no"} or \code{"yes"}. 
  }
  \item{coding}{
\code{"all"}, \code{"variable"}, \code{"noabsencesites"}, or \code{"nopresencesites"}. 
  }
  \item{parsmodel}{
\code{"no"} or \code{"yes"}.
  }
}
}
\references{
J.P. Huelsenbeck & Ronquist F. 2001. MrBayes: Bayesian inference of phylogenetic trees. \emph{Bioinformatics} \bold{17}: 754-755.

Ronquist F. & J.P. Huelsenbeck. 2003. MrBayes 3: Bayesian phylogenetic inference under mixed models. \emph{Biometrics} \bold{19}: 1572-1574.
	
MrBayes website: \url{https://mrbayes.sourceforge.net/}.
}
\author{Christoph Heibl}
%\note{}
\seealso{
\code{\link{mrbayes.prset}} to set prior distributions, \code{\link{mrbayes.mcmc}} to set parameters of the Markov chain Monte Carlo (MCMC), and \code{\link{mrbayes}} to run MrBayes locally or prepare input files for a computer cluster.
}
\examples{
## F81
mrbayes.lset(nst = 2)

## GTR + Gamma
mrbayes.lset(nst = 6, rates = "gamma")

## GTR + Gamma + I
mrbayes.lset(nst = 6, rates = "invgamma")
}