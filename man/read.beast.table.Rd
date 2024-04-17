\name{read.beast.table}
\Rdversion{1.1}
\alias{read.beast.table}

\title{
Extract node data from BEAST chronogram
}
\description{
This function reads a BEAST chronogram such as produced by TreeAnnotator and extracts time, rate, and support values for internal and external nodes. Nodes in the resulting data frame are ordered exactly like in the NEXUS file.
}
\usage{
read.beast.table(file, digits = 2)
}

\arguments{
  \item{file}{character string giving the input file, which must be a TreeAnnotator-generated chronogram in NEXUS format}
  
  \item{digits}{NULL or integer, if \code{!is.null(digits)} values are rounded to the given integer}
}

\value{
A matrix; each row corresponds to an internal node, the (ape!)number of which is given in the first column; the remaining columns list the node values extracted from the chronogram.
}

\author{Christoph Heibl}

\seealso{
\code{\link{read.beast}} to parse TreeAnnotator output, \code{\link{rbeauti}} to create XML input for BEAST. \code{HPDbars} for plotting highest posterior densities on phylogenies has been moved to package \bold{viper}.
}