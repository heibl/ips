\name{multistate}
\Rdversion{1.1}
\alias{multistate}
\alias{multistateMCMC}
\alias{multistateML}

\title{MULTISTATE}

\description{These functions provide wrappers to \code{BayesMultiState} in the BayesTraits package written by Mark Pagel and Andrew Meade.
}

\usage{
multistateML(phy, traits, model = "ARD", anc.states = TRUE, 
  path = "/Applications/BayesTraits", dir = NULL)

multistateMCMC(phy, traits, model = "ARD", anc.states = TRUE, 
    rd = 2, rjhp = NULL, fixNodes = NULL, it = 1e+05, bi = 10000,
    sa = 1000, path = "/Applications/BayesTraits", dir = NULL)
}

\arguments{
  \item{phy}{an object of class \code{phylo}.
}
  \item{traits}{a \code{data.frame} with two colums. The first column contains the taxon labels, the second column contains the character states.
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}

  \item{anc.states}{either \code{logical} or a \code{list}, the latter containing the tip labels of those internal nodes, for which the likelihood of ancestral character states should be estimated.
}
  \item{rd}{a real number, giving the \code{RateDev} parameter, i.e., the deviation of the normal distribution, that changes to the rates are drawn from. Should be set such that acceptance of the rate parameters is about 0.2.
}

  \item{rjhp}{a character string giving the details of priors and hyperpriors for the reversible jump MCMC (rjMCMC). If left \code{NULL}, a conventional MCMC is used. In order to use the rjMCMC, you must specify the distribution of the prior and the interval of the uniform hyperprior distribution that seeds it. For example, \code{exp 0 30} specifies an exponential distribution seeded from a uniform distribution on the interval 0 to 30, and \code{gamma 0 10 0 10} specifies a gamma prior with its mean and standard deviation seeded from uniform distributions on the interval 0 to 10.
}
  \item{fixNodes}{a list giving fixed character states of certain internal nodes. This argument corresponds to the \code{fossil} command in the MultiState manual.
}

  \item{it}{numeric, sets the number of iterations to run the MCMC for.
}
  \item{bi}{numeric, sets the number of iterations of the MCMC that will be discarded as burn-in.
}
  \item{sa}{numeric, sets the the sample period in the MCMC.
}
  \item{path}{a character string giving the path to executables in the BayesTraits package.
}
  \item{dir}{a character string giving a directory name where the input and output files will be stored. The directory will be created by \code{multistateML} and must not exist already. If \code{dir = NULL} (default) input and output is written to the working directory (thereby overwriting existing output).
}
}

\references{
The BayesTraits manual: \url{http://www.evolution.rdg.ac.uk/Files/BayesTraits-V1.0-Manual.pdf}.

Pagel, M., A. Meade, and D. Barker. 2004. Bayesian estiamation of ancestral character states on phylogenies. \emph{Syst. Biol.} \bold{53}: 673-684.

Pagel, M. and A. Meade. 2006. Bayesian analysis of correlated evolution of discrete characters by reversible-jump Markov chain Monte Carlo. \emph{Am. Nat.} \bold{167}: 808-825.
}

\author{Christoph Heibl}

\seealso{
	\code{\link{ace}}
}
