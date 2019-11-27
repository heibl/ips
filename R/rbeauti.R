## This code is part of the ips package
## © C. Heibl 2014 (last update 2016-12-05)

## to do: taxonsets
## to do: clocks

#' @title XML Input Files for BEAST
#' @description Prepare XML files for BEAST with R. BEAST uses an MCMC approach
#'   to estimate rooted phylogenies from molecular data (Drummond & Rambaut,
#'   2007).
#' @param ... One or more object(s) of class \code{\link{DNAbin}}.
#' @param file A connection, or a character string naming the file to write to.
#'   If left empty the XML tree will be printed to the screen (see Examples).
#' @param template \emph{Currently unused.}
#' @param taxonset A list containing one or more taxon sets.
# \item{monophyly}{A vector indicating monophyly constraints for the taxon sets declared with \code{taxonset}.}
# 
#  \item{tmrcaCons}{A list, containing the prior distribution(s) for age constraints of internal nodes (which must be grouped by \code{taxonset}).}
#  
#  \item{startingTree}{Either "random" or "upgma", or an object of class \code{"phylo"} to be used as a starting tree.}
#  
#  \item{specModel}{A character string indicating an evolutionary model to construct a prior distribution of node heights (tree prior). Currently implemented are the \bold{Yule Model} (\code{"yule"}), the \bold{Birth-Death Model} (\code{"birthDeath"}), and the \bold{Coalescent model with constant size} (choosen with any other string).}
#  
#  \item{clock}{A character string, either \code{"strict"} or \code{"lognormal"} to choose between the strict clock and the uncorrelated lognormal relaxed clock model. The exponential rates relaxed clock is currently not supported.}
#  
#  \item{ngen}{A character string, the number of generations to run the MCMC.}
#  
#  \item{samplefreq}{A character string, the intervals between sampling the MCMC.}
#  
#  \item{logSubs}{A logical, indicating if trees with branch lengths expressed as substitution should be logged.}
#  
#  \item{nodata}{A logical, indicating if BEAST should be run without data (see details).} 
#' @details \code{rbeauti} has been completely rewritten to work with
#'   \bold{BEAST 2}. Currently \code{rbeauti} offers few options, because the
#'   idea is not to create ready-to-use XML file. That can be done conveniently
#'   with \bold{BEAUti} (the BEAST package's genuine XML generator). Instead,
#'   \code{rbeauti} is intended to make the definition of large numbers of taxon
#'   sets easy. The creation of taxon sets can be done via R scripts and the
#'   resulting XML files can be further modified with BEAUti.
#' @references The BEAST 2 website:
#'   \url{http://beast.bio.ed.ac.uk/BEAST_v1.5.x_XML_Reference} Drummond, A.J. &
#'   A. Rambaut. 2007. BEAST: Bayesian evolutionary analysis by sampling trees.
#'   \emph{BMC Evolutionary Biology} \bold{7}: 240.
#' @seealso \code{\link{read.beast}}, \code{\link{read.beast.table}}
#' @examples 
#' data(ips.16S)
#' 
#' ## define taxon sets
#' spec <- rownames(ips.16S)
#' ingroup <- spec[grep("Ips|Orthomotomicus", spec)]
#' outgroup <- spec[grep("Pityogenes", spec)]
#' ts <- list(list(id = "ingroup", taxon = ingroup),
#'            list(id = "outgroup", taxon = outgroup))
#' 
#' ## print XML file to screen
#' rbeauti(ips.16S, taxonset = ts)
#' @import XML
#' @export

rbeauti  <- function(..., file, template = "standard", 
                     taxonset){
	
  ## handle partitions
  ## -----------------
  s <- list(...)
  if ( unique(sapply(s, class)) == "list" )
    s <- unlist(s, recursive = FALSE)
  
  ## assemble node(s) <data>
  ## -----------------------
  if ( is.null(names(s)) ){
    id <- paste("part", 1:length(data), sep = "")
    names(s) <- id
  } else {
    id <- names(s)
  }
  data <- assembleDataNode(s)
  
  ## assemble node <state>
  ## ---------------------
  state <- assembleStateNode(id)
  
  ## assemble node <init>
  ## ---------------------
  init <- lapply(id, assembleInitNode)
  
  ## assemble node <distribution>
  ## ----------------------------
  distribution <- assembleDistributionNode(id)
  
  ## assemble <operator> nodes
  ## -------------------------
  operators <- assembleOperators(id)
  
  ## assemble node <run> 
  ## ---------------------
  run <- xmlNode("run", 
                 attrs = c(chainLength = "10000000", 
                                  id = "mcmc",
                                  spec = "MCMC"))
  run <- addChildren(run, kids = list(state))
  run <- addChildren(run, kids = init)
  run <- addChildren(run, kids = list(distribution))
  run <- addChildren(run, kids = operators)
  run <- addChildren(run, kids = assembleLoggers(id))
  
  ## assemble nodes <map> 
  ## --------------------
  m <- mm <- c("Beta", "Exponential", "InverseGamma", "LogNormalDistributionModel", "Gamma", 
            "Uniform", "Prior", "LaplaceDistribution", "OneOnX", "Normal")
  m[m == "Prior"] <- "prior"; m[m == "LogNormalDistributionModel"] <- "LogNormal"
  mm <- paste("beast.math.distributions", mm, sep = ".")
  m <- cbind(name = m, mm)
  maps <- apply(m, 1, function(x) xmlNode("map", x[2], 
                                     attrs = x[1]))
  
  
  ## assemble node <beast> 
  ## ---------------------
  namespace <- c("beast.core", "beast.evolution.alignment", "beast.evolution.tree.coalescent", 
                 "beast.core.util", "beast.evolution.nuc", "beast.evolution.operators", 
                 "beast.evolution.sitemodel", "beast.evolution.substitutionmodel", "beast.evolution.likelihood")
  namespace <- paste(namespace, collapse = ":")
  beast <- xmlNode("beast", attrs = c(beautitemplate = "Standard",
                                      beautistatus = "",
                                      namespace = namespace,
                                      version = "2.0"))
  beast <- addChildren(beast, kids = data)
  beast <- addChildren(beast, kids = maps)
  beast <- addChildren(beast, kids = list(run))
  ## convert from class XMLNode to XMLInternalDocument
  ## -------------------------------------------------
#   beast <- saveXML(beast, encoding = "UTF-8",
#                    prefix = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
#   beast <- xmlInternalTreeParse(beast, asText = TRUE)
  
  if ( missing(file) ) return(beast)
  else {
    if ( length(grep("[.]xml$", file)) == 0 ) 
      file <- paste(file, "xml", sep = ".")
    saveXML(beast, file = file, encoding = "UTF-8",
            prefix = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
  }
}