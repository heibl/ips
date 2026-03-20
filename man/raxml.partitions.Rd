\name{raxml.partitions}
\alias{raxml.partitions}
\docType{data}
\title{Partition scheme for RAxML
}
\description{
Given a set of DNA sequence alignments, \code{raxml.partitions} creates a data frame with partition bounderies that can be input into \code{\link{raxml}}.
}
\usage{raxml.partitions(...)}
\arguments{
\item{...}{Two or more DNA sequence alignments of class \code{\link[ape]{DNAbin}}.}
}
\details{
For \code{raxml.partitions} to make sense, the DNA sequence alignments must be given exactly in the same order in which they are concatenated into a supermatrix (see Examples section). Without any testing, the type of sequences is supposed to be DNA.
}
\value{
A data frame with four columns (\code{type}, \code{locus}, \code{begin}, and \code{end}) and number of rows corresponding to the number of partitions.
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\seealso{
\code{\link[ape]{cbind.DNAbin}} to concatenate multiple alignments; \code{\link{raxml}} for an interface to RAxML.
}
\examples{
## bark beetle sequences
data(ips.cox1)
data(ips.16S)
data(ips.28S)

## Note the same order of individual
## alignments in both functions:
## ----------------------------
raxml.partitions(cox1 = ips.cox1, 
                 r16S = ips.16S, 
                 r28S = ips.28S)

cbind(ips.cox1, ips.16S, ips.28S,
  fill.with.gaps = TRUE)
}
