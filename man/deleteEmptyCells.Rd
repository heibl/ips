\name{deleteEmptyCells}
\alias{deleteEmptyCells}

\title{Delete Spurious Rows and Columns from DNA Alignments}

\description{
After subsetting (see e.g. \code{\link[ape]{DNAbin}}), DNA sequence alignments can contain rows and columns that consist entirely of missing and/or ambiguous character states. \code{deleteEmptyCells} will delete all such rows (taxa) and columns (characters) from a DNA sequence alignment.
}

\usage{
deleteEmptyCells(DNAbin, nset = c("-", "n", "?"), 
  quiet = FALSE)
}

\arguments{
  \item{DNAbin}{an object of class \code{\link[ape]{DNAbin}}.
  }
  \item{nset}{vector of mode character; rows or columns that consist \bold{only} of the characters given in \code{nset} will be deleted from the alignment.
  }
  \item{quiet}{logical: if set to \code{TRUE}, screen output will be suppressed}
}

\details{
For faster execution, \code{deleteEmptyCells} handles sequences in \bold{ape}'s bit-level coding scheme. As a consequence, \code{nset} cannot be extended by the user and the use of IUPAC ambiguity symbols other than \code{'N'} is currently not possible.
}

\value{
  an object of class \code{\link[ape]{DNAbin}}.
}

\author{Christoph Heibl}

\seealso{
  \code{\link{trimEnds}}, \code{\link{deleteGaps}}
}
\examples{
# COX1 sequences of bark beetles
data(ips.cox1)

# introduce completely ambiguous rows and columns
x <- as.character(ips.cox1[1:6, 1:60])
x[3, ] <- rep("n", 60)
x[, 20:24] <- rep("-", 6)
x <- as.DNAbin(x)
image(x)

# delete those rows and columns
x <- deleteEmptyCells(x)
image(x)
}
