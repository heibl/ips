\name{write.partitioned.nex}
\Rdversion{1.1}
\alias{write.partitioned.nex}

\title{Write partitioned and commented NEXUS files}

\description{
This function takes one or more objects of class \code{"DNAbin"} and writes them to an interleaved NEXUS file.
}

\usage{
write.partitioned.nex(..., file, labels)
}

\arguments{
  \item{\dots}{One or more objects of class \code{"DNAbin"}}
  \item{file}{character string giving a file name
}
  \item{labels}{vector of mode \code{"character"; used to label partitions; its length must match the number of gene partitions}
}
}
 
\author{Christoph Heibl}

\seealso{\code{\link{write.nex}}, \code{\link{write.phy}}, \code{\link{write.fas}}}

\examples{
# load sequence data and split in three parts
# -------------------------------------------
data(woodmouse)
w1 <- woodmouse[, 1:250]
w2 <- woodmouse[, 251:750]
w3 <- woodmouse[, 750:965]

# write nexus file with taxon block
# ---------------------------------
markers <- c("gene1", "gene2", "gene3")
write.partitioned.nex(w1, w2, w3, labels = markers)
}
