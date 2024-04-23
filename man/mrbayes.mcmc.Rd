\name{mrbayes.mcmc}
\alias{mrbayes.mcmc}
\title{MCMC Settings for MrBayes}
\description{
Set Markov chain Monte Carlo (MCMC) parameters for \code{\link{mrbayes}}.
}
\usage{
mrbayes.mcmc(...)
}
\arguments{
  \item{\dots}{
arguments in \code{tag = value} form, or a list of tagged values. The tags must come from the names of MCMC parameters described in the ‘MCMC Parameters’ section.
  }
}
%\details{}
\value{
a list containing a subset (including the empty and the full set) of model parameters.
}
\section{MCMC Parameters}{
\describe{
 \item{ngen}{ 
"NUMERIC"
  }
  \item{nruns}{ 
"NUMERIC"
  }
  \item{nchains}{ 
"NUMERIC"
  }
  \item{temp}{ 
"NUMERIC"
  }
%    # reweight        <number>,<number>     0.00 v 0.00 ^,
  \item{swapfreq}{ 
"NUMERIC"
  }
  \item{nswaps}{ 
"NUMERIC"
  }
  \item{samplefreq}{ 
"NUMERIC"
  }
  \item{printfreq}{ 
"NUMERIC"
  }
  \item{printall}{ 
\code{"yes"} or \code{"no"}
  }
  \item{printmax}{ 
"NUMERIC"
  }
  \item{mcmcdiagn}{ 
\code{"yes"} or \code{"no"}
  }
  \item{diagnfreq}{ 
"NUMERIC"
}
  \item{diagnstat}{ 
\code{"avgstddev"} or \code{"maxstddev"}
  }
  \item{minpartfreq}{ 
"NUMERIC"
  }
  \item{allchains}{ 
\code{"yes"} or \code{"no"}
  }
  \item{allcomps}{ 
\code{"yes"} or \code{"no"}
  }
  \item{relburnin}{ 
\code{"yes"} or \code{"no"}
  }
  \item{burnin}{ 
"NUMERIC"
  }
  \item{burninfrac}{ 
"NUMERIC"
  }
  \item{stoprule}{ 
\code{"yes"} or \code{"no"}
  }
  \item{stopval}{ 
"NUMERIC"
  }
  \item{savetrees}{ 
\code{"yes"} or \code{"no"}
  }
  \item{checkpoint}{ 
\code{"yes"} or \code{"no"}
  }
  \item{checkfreq}{ 
"NUMERIC"
  }
%   # filename = filename,
  \item{startparams}{ 
\code{"current"} or \code{"reset"}
  }
  \item{starttree}{ 
\code{"current"}, \code{"random"}, or \code{"parsimony"}
  }
  \item{nperts}{ 
"NUMERIC"
  }
  \item{data}{ 
\code{"yes"} or \code{"no"}
  }
  \item{ordertaxa}{ 
\code{"yes"} or \code{"no"}
  }
  \item{append}{ 
\code{"yes"} or \code{"no"}
  }
  \item{autotune}{
\code{"yes"} or \code{"no"}
  }
  \item{tunefreq}{ 
"NUMERIC"
  }
}
}
\references{
J.P. Huelsenbeck & Ronquist F. 2001. MrBayes: Bayesian inference of phylogenetic trees. \emph{Bioinformatics} \bold{17}: 754-755.

Ronquist F. & J.P. Huelsenbeck. 2003. MrBayes 3: Bayesian phylogenetic inference under mixed models. \emph{Biometrics} \bold{19}: 1572-1574.
	
MrBayes website: \url{https://mrbayes.sourceforge.net/}.
}
\author{Christoph Heibl}
\note{The parameters \bold{reweight} and \bold{filename} cannot be set via \code{mrbayes.mcmc}.}
\seealso{
\code{\link{mrbayes.lset}} to set model parameters, \code{\link{mrbayes.prset}} to set prior distributions, and \code{\link{mrbayes}} to run MrBayes locally or prepare input files for a computer cluster.
}
\examples{
mrbayes.mcmc()
}
