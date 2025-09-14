## This code is part of the ips package
## © C. Heibl 2014 (last update 2025-09-14)

## to do: taxonsets
## to do: clocks

#' @importFrom XML addChildren saveXML
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