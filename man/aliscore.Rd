\name{aliscore}
\alias{aliscore}
\title{Masking of Sequence Alignments with ALISCORE}

\description{
This function provides a wrapper to \bold{Aliscore}, which can be used remove problematic regions of a sequence alignment.
}

\usage{
aliscore(x, gaps = "5state", w = 6, r, t, l, s, o, 
        path = "/Applications/Aliscore_v.2.0")
}

\arguments{
  \item{x}{DNA sequences of class \code{DNAbin}.}
  
  \item{gaps}{A vector of mode \code{"character"} indicating how gaps shall be treated: as \code{"5state"} or as \code{"ambiguous"}.}
  
  \item{w}{An integer giving the size of the sliding window.}
  
  \item{r}{An integer giving the number of random pairwise sequence comparisons; defaults to 4*N.}
  
  \item{t}{\emph{Not yet implemented}.}
  
  \item{l}{\emph{Not yet implemented}.}
  
  \item{s}{\emph{Not yet implemented}.}
  
  \item{o}{A vector of mode \code{"character"} containing outgroup taxon names.}

  \item{path}{A character string, giving the path to the Aliscore script.}
}

\value{
    \code{matrix} of class \code{"DNAbin"}
}

\references{
Misof, B. and K. Misof. 2009. A Monte Carlo approach successfully identifies randomness in multiple sequence alignments: a more objective means of data exclusion. \emph{Syst. Biol.} \bold{58}: 21--34.

Kueck, P., K. Meusemann, J. Dambach, B. Thormann, B.M. von Reumont, J.W. Waegele and B. Misof. 2010. Parametric and non-parametric masking of randomness in sequence alignments can be improved and leads to better resolved trees. \emph{Frontiers in Zoology} \bold{7}: 10.
}

\author{
Christoph Heibl
}

\note{
This function was developed with ALISCORE version 2. 
}

\seealso{
\code{\link{mafft}} and \code{\link{prank}} for sequence alignment; \code{\link{gblocks}} for another alignment masking algorithm.
}
%\example{
%# aliscore(x)
%}
