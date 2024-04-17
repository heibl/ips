\name{trimEnds}
\alias{trimEnds}
\alias{fillEndsWithN}
\title{Trim Alignment Ends}

\description{
This function trims both ends of a DNA sequence alignment to the first and last alignment positions that contain a minimum number of non-ambiguous characters (\code{a}, \code{c}, \code{g}, \code{t}). In addition, all gap characters (\code{"-"}) beyond the first and last non-ambiguous characters of each sequence are replaced by the  character \code{"n"}.
}

\usage{
trimEnds(x, min.n.seq = 4)
}

\arguments{
  \item{x}{An object of class \code{DNAbin}.}
  \item{min.n.seq}{A \code{numeric} giving the required minimum number of sequences having an non-ambiguous base character (a, c, g, t) in the first and last position of the alignment; defaults to \code{4}, which is the minimum number of sequences needed to produce a non-trivial unrooted topology.}
}

\value{
  An object of class \code{DNAbin}.
}

\author{Christoph Heibl}

\seealso{
\code{\link{deleteEmptyCells}}, \code{\link{deleteGaps}}
}

\examples{
# simple example alignment:
# -------------------------
x <- structure(list(nb = 5, 
	seq = c("acaaggtaca", "-caaggtac-", "acaaggtaca", "aca--gtaca", "-ccaggta--"), 
	nam = LETTERS[1:5]), .Names = c("nb", "seq", "nam"), 
	class = "alignment")

# convert to DNAbin:
# ------------------
x <- as.DNAbin(x)

# fill missing nucleotides:
# -------------------------
x <- trimEnds(x)

# show results:
# -------------
as.character(x[2, ])
}
