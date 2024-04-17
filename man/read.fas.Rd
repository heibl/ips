\name{read.fas}
\alias{read.fas}
\alias{read.nex}
\alias{read.phy}

\title{Read DNA Sequences}

\description{
This functions parse DNA sequences in FASTA, PHILIP, and NEXUS formatted files.
}

\usage{
read.fas(x, text)

read.nex(x)

read.phy(x)
}

\arguments{
  \item{x}{a character string, giving the file name.}
  \item{text}{a character string in FASTA format.}
}

\value{
  An object of class \code{DNAbin} and of mode \code{matrix} if all sequences are of same length or \code{matrix} if they are not.
}

\references{
Maddison, D.R., D.L. Swofford, and W.P. Maddison. 1997. NEXUS: an extensible file format for systematic information. \emph{Syst. Biol.} \bold{46}: 590-621.

}

\author{Christoph Heibl}

\seealso{\code{\link{mafft}} and \code{\link{prank}} for sequence alignment, \code{\link{gblocks}} and \code{\link{aliscore}} for quality check and cleaning of sequence alignments, \code{\link{cbind.DNAbin}} for concatenation of sequence alignments.}

\examples{
## bark beetle COX1 sequences
## --------------------------
data(ips.cox1)

## create temporary file names
## ---------------------------
format <- c(".fas", ".phy", ".nex")
fn <- sapply(format, tempfile, 
  pattern = "ips", tmpdir = tempdir())

## write sequences files
## ---------------------
write.fas(ips.cox1, fn[".fas"])
write.phy(ips.cox1, fn[".phy"])
write.nex(ips.cox1, fn[".nex"])

## read sequence files
## -------------------
fas <- read.fas(fn[".fas"])
phy <- read.phy(fn[".phy"])
nex <- read.nex(fn[".nex"])

## remove sequence files
## ---------------------
unlink(fn)
}
