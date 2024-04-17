\name{mrbayes}
\Rdversion{1.1}
\alias{mrbayes}
\title{Bayesian MCMC Tree Search with MrBayes}
\description{
Provides a wrapper for Bayesian phylogenetic tree search through MrBayes (Huelsenbeck & Ronquist, 2001; Ronquist & Huelsenbeck, 2003).
}
\usage{
mrbayes(x, file = "", lset, prset, mcmc, unlink, constraint,
<<<<<<< HEAD
        burnin = 10, contype = "allcompat", run = FALSE)
}
\arguments{
  \item{x}{
an object of class \code{\link{DNAbin}} in the case of \code{mrbayes} or a matrix of mode character in the case of \code{mrbayes.mixed}.
  }
  \item{file}{
a character string, giving the name of the MrBayes input file.
  }
  \item{lset}{
a list as returned by \code{\link{mrbayes.lset}} containing the parameter settings of the model of molecular evolution.  
  }
  \item{prset}{
 a list as returned by \code{\link{mrbayes.prset}} containing the parameter setting for the prior distirbutions. 
  }
  \item{mcmc}{
a list as returned by \code{\link{mrbayes.mcmc}} containing the parameter setting for the Marcov chain Monte Carlo (MCMC).  
=======
        burnin = 10, contype = "allcompat", exec, run = FALSE)
}
\arguments{
  \item{x}{
An object of class \code{\link{DNAbin}} in the case of \code{mrbayes}.
  }
  \item{file}{
A character string, giving the name of the MrBayes input file.
  }
  \item{lset}{
A list as returned by \code{\link{mrbayes.lset}} containing the parameter settings of the model of molecular evolution.  
  }
  \item{prset}{
 A list as returned by \code{\link{mrbayes.prset}} containing the parameter setting for the prior distributions. 
  }
  \item{mcmc}{
A list as returned by \code{\link{mrbayes.mcmc}} containing the parameter setting for the Markov chain Monte Carlo (MCMC).  
>>>>>>> 647ba4c9ae1104b649cacc901c527de0207b7775
  }
  \item{unlink}{
  
  }
  \item{constraint}{
  
  }
  \item{burnin}{
<<<<<<< HEAD
an integer; the number of samples from the MCMC to be discarded prior to further analysis.
  }
  \item{contype}{
a character string; the type of consensus tree calculated from the posterior distribution of trees: either \code{"halfcompat"} (majority-rule consensus tree) or \code{"allcombat"} (strict consensus tree).
  }
  \item{run}{
logical; \code{run = FALSE} will only print the NEXUS file, \code{run = TRUE} will also start the MCMC runs, if the \code{path} argument is correctly specified.
=======
An integer; the number of samples from the MCMC to be discarded prior to further analysis.
  }
  \item{contype}{
A character string; the type of consensus tree calculated from the posterior distribution of trees: either \code{"halfcompat"} (majority-rule consensus tree) or \code{"allcombat"} (strict consensus tree).
  }
  \item{exec}{
  A character string giving the full path of the MrBayes program.
  }
  \item{run}{
Logical; \code{run = FALSE} will only print the NEXUS file, \code{run = TRUE} will also start the MCMC runs, if \code{exec} is correctly specified.
>>>>>>> 647ba4c9ae1104b649cacc901c527de0207b7775
  }
}
\details{ 	
\code{mrbayes} was last updated and tested with MrBayes \bold{v3.2.2} under R 3.1.0 on a x86_64-apple-darwin10.8.0 (64-bit) platform. It is intended to offer a simply parameterized building block for larger scripts.	
}

\value{
	None; a NEXUS file with MrBayes block is written to a file and, if \code{run = TRUE}, the MCMC runs in MrBayes are started.
}

\references{
J. P. Huelsenbeck & Ronquist F. 2001. MrBayes: Bayesian inference of phylogenetic trees. \emph{Bioinformatics} \bold{17}: 754-755.

Ronquist F. & J. P. Huelsenbeck. 2003. MrBayes 3: Bayesian phylogenetic inference under mixed models. \emph{Biometrics} \bold{19}: 1572-1574.
	
MrBayes website: \url{http://mrbayes.sourceforge.net/}.
}
\author{Christoph Heibl}

\seealso{
\code{\link{mafft}} and \code{\link{prank}} for sequence alignment; \code{\link{raxml} for maximum likelihood tree search.}
}

\examples{	
data(ips.cox1)
x <- ips.cox1[, 100:140] # tiny alignment
mrbayes(x, file = "", mcmc = mrbayes.mcmc(ngen = 100), run = FALSE)

\dontrun{

library(phangorn)
tree <- rtree(10)
Y1 <- simSeq(tree, l = 20)
Y2 <- simSeq(tree, l = 20, type = "USER", levels=c("0", "1"))
Y <- cbind(as.character(Y1), as.character(Y2))
<<<<<<< HEAD
mrbayes.mixed(Y, file = "", ngen = 100, run = FALSE)
=======
mrbayes(Y, file = "", run = FALSE)
>>>>>>> 647ba4c9ae1104b649cacc901c527de0207b7775
}

}