% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rc.R
\name{rc}
\alias{rc}
\title{Reverse-Complement of DNA sequences}
\usage{
rc(seqs, i, complement = TRUE, reverse = TRUE, delete.gaps = FALSE)
}
\arguments{
\item{seqs}{An object of class \code{DNAbin}.}

\item{i}{Logical or numeric index to indicate a subset of sequences to
manipulate. If not given the entire set of sequences will be manipulated.}

\item{complement}{Logical, indicating if sequences will be turned into their
complement.}

\item{reverse}{Logical, indication if sequences will be reversed.}

\item{delete.gaps}{Logical, indicating if gap character will be removed prior
to sequence manipulation.}
}
\value{
An object of the same class and dimension as \code{seqs}.
}
\description{
Reverse, complement or reverse-complement of DNA sequences.
}
\examples{
## A minimal sequence alignment:
x <- list(
  seqA = c("a", "a", "c", "c", "g", "t"),
  seqB = c("n", "-", "r", "y", "g", "t"))
x <- as.DNAbin(x)

## Three choices of manipulation:
as.character(x)
as.character(rc(x))                      ## reverse-complement
as.character(rc(x, complement = FALSE))  ## only reverse
as.character(rc(x, reverse = FALSE))     ## only complement

## You can remove gaps:
as.character(rc(x, delete.gaps = TRUE))  ## gaps/indels removed

}
