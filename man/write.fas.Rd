\name{write.fas}
\alias{write.fas}
\alias{write.phy}
\alias{write.nex}

\title{Write DNA Sequences to File}

\description{
Write DNA sequences and morphological data to FASTA, PHYLIP, or NEXUS formatted files.
}

\usage{
	
write.fas(x, file, block.width = FALSE, 
    truncate = FALSE, append = FALSE)

write.phy(x, file, block.width = FALSE, 
    strict = FALSE)

write.nex(x, file, block.width = 60, 
    taxblock = FALSE)
}

\arguments{
  \item{x}{
an object of class \code{\link{DNAbin}} (usually as matrix, but \code{write.fas} also accepts lists) or a list of objects of class \code{DNAbin} (only \code{write.nex}) or a data frame containing standard (morphological, etc.) data (only \code{write.nex}).
  }
  \item{file}{
a character string giving the filename; a special case is \code{file = ""}, which causes the file content to be written on the standard output connection (i.e. the console). If \code{file} is left unspecified (default), the file content is returned as a vector of mode \code{"character"} and can be used as a building block for more complex data files.
  }
  
  \item{block.width}{
an integer, giving the number of characters per line.
  }
  \item{truncate}{
truncation of taxon names to the number of characters given as a integer, otherwise (default) taxon names will not be changed.
  }
  \item{append}{
logical, if \code{TRUE} the sequences will be appended to \code{file} (if it exists).
  }
  \item{strict}{
logical, if \code{TRUE} the names of the sequences will be truncated to 10 strings.
  }
  \item{taxblock}{
logical, if \code{TRUE}, a tax block will be added to the NEXUS file.
  }
}
\details{
\code{write.nex} can handle multiple DNA sequence alignments, which are handed over as a list of objects of class \code{DNAbin}. Correct matching of the rows in the alignments is cared for automatically, hence the individual alignments can contain different numbers of samples and samples need not be in the same order.
}
\value{
None, except when called with \code{file} left unspecified, which causes the file content to be returned as a vector of mode \code{"character"}. This is particularly useful for constructing special types of input files, e.g. for MrBayes (\code{\link{mrbayes}}).
}

\references{
Maddison, D.R., D.L. Swofford, and W.P. Maddison. 1997. NEXUS: an extensible file format for systematic information. \emph{Syst. Biol.} \bold{46}: 590-621.}
 
\author{Christoph Heibl}

\seealso{
	\code{\link{read.fas}}, \code{\link{read.phy}}, and \code{\link{read.nex}} for reading of DNA sequence files.
}

\examples{
data(ips.cox1)
data(ips.28S)

## Examples for FASTA files
## ------------------------
write.fas(ips.cox1[1:5, 1:120], block.width = 60)

## Examples for PHYLIP files
## ------------------------
write.phy(ips.cox1[1:5, 1:20], block.width = 40)

## Examples for NEXUS files
## ------------------------
x <- list(cox1 = ips.cox1[1:5, 1:10], 
          rna28S = ips.28S[1:5, 1:30])
write.nex(x, block.width = 20)

# Truncation of taxonnames:
# -------------------------
rownames(ips.cox1)[1] <- "AVeeeeeeeeeeeeeeeeeryLongName"
write.fas(ips.cox1, truncate = 10)

# If truncation leads to identical taxonnames,
# a warning will be issued:
# -------------------------
rownames(ips.cox1)[1:2] <- "AVeeeeeeeeeeeeeeeeeryLongName"
write.fas(ips.cox1, truncate = 10)
}