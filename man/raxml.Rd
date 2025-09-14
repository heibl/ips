\name{raxml}
\alias{raxml}

\title{Maximum Likelihood Tree Estimation with RAxML}

\description{
This function calls \bold{RAxML} (see Reference section) for the maximum likelihood estimation of tree topology and/or branch lengths for a given DNA sequence alignment, rapid and conventional non-parametric bootstrapping, mapping splits onto individual topologies, and a lot more. See the RAxML manual for details, especially if you are a new user of RAxML. 
}

\usage{
raxml(DNAbin, m, f, N, p, b, x, k, 
      partitions, outgroup, backbone = NULL, 
      file = "fromR", exec, threads)
}

\arguments{
  \item{DNAbin}{a matrix of DNA sequences of class \code{\link[ape]{DNAbin}}.}
  
  \item{m}{a vector of mode \code{"character"} choosing a model of substitution; currently only GTR models available}
  
  \item{f}{a vector of mode \code{"character"} choosing from the RAxML algorithm analogous to the \code{-f} flag (see Detail section and RAxML manual).}
  
  \item{N}{either of mode \code{"integer"} or \code{"character"}. Integers give the number of independent searches on different starting tree or replicates in bootstrapping. Alternatively, one of four bootstopping criteria can be chosen: \code{"autoFC"}, \code{"autoMR"}, \code{"autoMRE"}, or \code{"autoMRE_IGN"}.
  }
  
  \item{p}{integer setting a random seed for parsimony starting trees.}
  
  \item{b}{integer setting a random seed for bootstrapping.}
  
  \item{x}{integer setting a random seed for rapid bootstrapping.}
  
  \item{k}{logical, if \code{TRUE} the branch lengths of bootstrapped trees are recorded}
  
  \item{partitions}{a \code{data frame} specifying partitions of the alignment.}
  
  \item{outgroup}{a vector of mode \code{"character"} containing the names of the outgroup taxa.
  }
  
  \item{backbone}{a \code{phylo} object handed as a backbone tree to RAxML.}
    
  \item{file}{a vector of mode \code{"character"} giving a name to the output files.}
  
  \item{exec}{a vector of mode \code{"character"} giving the path to the directory containing the RAxML executable. The default value will work on Mac OS X if the folder containing the executable is renamed to \code{"RAxML-8.0.3"}.}
  
  \item{threads}{integer giving the number of parallel threads to use (PTHREADS only).}
  
}

\details{
There are some limitations of this wrapper compared to RAxML run directly from the command line. 
\enumerate{
\item Only DNA is allowed as data type.
\item Option \code{f} can only take a limited number of values (\code{d}, \code{a}).
} % close enumerate

RAxML needs the specification of random seeds for parsimony estimation of starting trees and for bootstrap resampling. The corresponding argument names in \code{raxml} are identical to the flags used by RAxML (\code{-p}, \code{-b}, and \code{-x}). If you choose not to give any values, \code{raxml} will generate a (different) value for each required random seed every time it is called. Be aware that \code{\link{set.seed}} will work only for \code{p}, but not for \code{b} or \code{x}.
} % close details

\value{
  A list with a variable number of elements, depending on the analysis chosen: 
  
  \tabular{ll}{
  \code{"info"} \tab RAxML log file as character string\cr
  \code{"bestTree"} \tab MLE of tree\cr
  \code{"bipartitions"} \tab MLE of tree annotated with bootstrap proportions\cr
  \code{"bootstrap"} \tab bootstrapped trees\cr
  }
  
}
% REFERENCES
\references{
(in chronological order)

Stamatakis, A., T. Ludwig and H. Meier. 2004. RAxML-III: A fast program for maximum likelihood-based inference of large phylogenetic trees. \emph{Bioinformatics} \bold{1}: 1--8.

Stamatakis, A. 2006. RAxML-VI-HPC: 
Maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models. \emph{Bioinformatics} \bold{22}: 2688--2690.

Stamatakis, A., P. Hoover, and J. Rougemont. 2008. A rapid bootstrap algorithm for the RAxML web-servers. \emph{Syst. Biol.} \bold{75}: 758--771.

Pattengale, N. D., M. Alipour, O. R. P. Bininda-Emonds, B. M. E. Moret, and A. Stamatakis. 2010. How many bootstrap replicates are necessary? \emph{Journal of Computational Biology} \bold{17}: 337-354.

Stamatakis, A. 2014. RAxML Version 8: A tool for phylogenetic analysis and post-analysis of large phylogenies. \emph{Bioinformatics} Advance Access.
}
% AUTHOR
\author{Christoph Heibl}
% NOTE
\note{
RAxML is a C program and the source code is not contained in this package. This means that in order to run this function you will need to install RAxML yourself. See the 'Software' tab on \url{http://www.exelixis-lab.org/} for the most recent documentation and source code of RAxML. Depending on where you chose to install RAxML, you need to adjust the \code{exec} argument.

\code{raxml} was last tested and running fine on Mac OS X with RAxML 8.0.29. Please be aware that calling third-party software from within R is a platform-specific process and I cannot guarantee that \code{raxml} will behave properly on any system. 
}
% SEEALSO
\seealso{
\code{\link{raxml.partitions}} to store partitioning information in a data frame suitable for input as \code{partitions} argument in \code{raxml}.
}
% EXAMPLES
\examples{
## bark beetle sequences
data(ips.cox1)
data(ips.16S)
data(ips.28S)

ips <- cbind(ips.cox1, ips.16S, ips.28S,
  fill.with.gaps = TRUE)

exec <- NULL # replace by your RAxML path

if ( !is.null(exec) ){

## normal tree search with GTRCAT and GTRGAMMA
tr <- raxml(woodmouse, f = "d", N = 2, p = 1234, 
            exec = exec) # -1743.528461
tr <- raxml(woodmouse, m = "GTRGAMMA", f = "d", N = 2, p = 1234, 
            exec = exec)

## rapid bootstrap
tr <- raxml(woodmouse, m = "GTRGAMMA",
            f = "a", N = 10, p = 1234, x = 1234, 
            exec = exec)

## rapid bootstrap with automatic halt
tr <- raxml(woodmouse, m = "GTRGAMMA", 
            f = "a", N = "autoMRE", p = 1234, x = 1234, 
            exec = exec)
}

}