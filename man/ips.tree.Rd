\name{ips.tree}
\alias{ips.tree}
\docType{data}
\title{Ips Phylogeny}
\description{Phylogenetic tree of bark beetles (genus \emph{Ips}).
}
\usage{data(ips.tree)}
\format{
  The format is:
List of 5
 $ edge       : int [1:72, 1:2] 38 39 39 40 41 42 42 43 44 45 ...
 $ Nnode      : int 36
 $ tip.label  : chr [1:37] "Ips_acuminatus" "Ips_duplicatus" "Ips_integer" "Ips_plastographus" ...
 $ edge.length: num [1:72] 0.2806 0.0727 0.0295 0.0097 0.021 ...
 $ node.label : chr [1:36] "" "100" "21" "12" ...
 - attr(*, "class")= chr "phylo"
 - attr(*, "order")= chr "cladewise"
}
%\details{}
%\source{}
%\references{}
\examples{
data(ips.tree)
plot(ips.tree)
}
\keyword{datasets}
