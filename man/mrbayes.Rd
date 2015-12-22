\name{mrbayes}
\Rdversion{1.1}
\alias{mrbayes}
\alias{mrbayes.mixed}

\title{Bayesian MCMC Tree Search with MrBayes}

\description{
This function is a wrapper for Bayesian phylogenetic tree search through MrBayes (Ronquist & Huelsenbeck, 2003) with either DNA (\code{mrbayes}) or morphological (\code{mrbayes.mixed}) data.
}

\usage{
mrbayes(x, file = "", nst = 6, rates = "invgamma", ngammacat = 4, 
    nruns = 2, ngen = 1e+06, printfreq = 100, samplefreq = 10, 
    nchains = 4, savebrlens = "yes", temp = 0.2, burnin = 10, 
    contype = "allcompat", run = FALSE)


mrbayes.mixed(x, file, nst = 6, rates = "invgamma", ngammacat = 4, nruns = 2, 
	ngen = 1e+06, printfreq = 100, samplefreq = 10, nchains = 4,
	savebrlens = "yes", temp = 0.2, burnin = 10, contype = "allcompat", run = TRUE)
}

\arguments{
  \item{x}{An object of class \code{\link{DNAbin}} in the case of \code{mrbayes} or a matrix of mode character in the case of \code{mrbayes.mixed}.}
  
  \item{file}{A character string, giving the name of the MrBayes input file.}
  
  \item{nst}{An integer giving the number of rates in the model of sequence evolution.}
  
  \item{rates}{A character string; allowed are \code{"equal"}, \code{"gamma"}, \code{"propinv"}, \code{"invgamma"}, and \code{"adgamma"}; the default is \code{"equal"}.}
  
  \item{ngammacat}{An integer; the number rate categories for the discretized Gamma distribution; the default is \code{4}.}
  
  \item{nruns}{An integer; the number of runs.}
  
  \item{ngen}{An integer; the number of states of the MCMC.}
  
  \item{printfreq}{An integer; the interval between states of the MCMC to be printed on the screen}
  
  \item{samplefreq}{An integer; the interval between states of the MCMC to be \bold{sampled}.
}

  \item{nchains}{An integer; number of Metropolis coupled MCMCs in each run.}
  
  \item{savebrlens}{Logical; shall branch lengths be saved.}
  
  \item{temp}{
%%     ~~Describe \code{temp} here~~
}

  \item{burnin}{An integer; the number of samples from the MCMC to be discarded prior to further analysis.}
  
  \item{contype}{A character string; the type of consensus tree calculated from the posterior distribution of trees: either \code{"halfcompat"} (majority-rule consensus tree) or \code{"allcombat"} (strict consensus tree).
}
  
  \item{run}{Logical; \code{run = FALSE} will only print the NEXUS file, \code{run = TRUE} will also start the MCMC runs, if the \code{path} argument is correctly specified.}
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
mrbayes(x, file = "", ngen = 100, run = FALSE)

%\dontrun{

library(phangorn)
tree <- rtree(10)
Y1 <- simSeq(tree, l = 20)
Y2 <- simSeq(tree, l = 20, type = "USER", levels=c("0", "1"))
Y <- cbind(as.character(Y1), as.character(Y2))
mrbayes.mixed(Y, file = "", ngen = 100, run = FALSE)
%}

}