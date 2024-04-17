## This code is part of the ips package
## © C. Heibl 2014 (last update 2020-03-31)

## to do: taxonsets
## to do: substitution models

#' @title XML Input Files for BEAST
#' @description Prepare XML files for BEAST with R. BEAST uses an MCMC approach
#'   to estimate rooted phylogenies from molecular data (Drummond & Rambaut,
#'   2007).
#' @param ... One or more object(s) of class \code{\link{DNAbin}}.
#' @param file A connection, or a character string naming the file to write to.
#'   If left empty the XML tree will be printed to the screen (see Examples).
#' @param template \emph{Currently unused.}
#' @param link.clocks Logical, indicating if clock models should be linked over partitions.
#' @param link.trees Logical, indicating if tree models should be linked over partitions.
#' @param subst.model A character string defining a substituion model, either
#'   \code{"JC69"}, \code{"HKY"},\code{"TN93"}, or \code{"GTR"}.
#' @param clock A character string defining a clock model, either \code{"Strict
#'   Clock"}, \code{"Relaxed Clock Exponential"},\code{"Relaxed Clock Log
#'   Normal"}, or \code{"Random Local Clock"}.
#' @param tree A character string defining a tree model.
#' @param taxonset A list containing one or more taxon sets.
#' @param chain.length Integer, the number of generations to run the MCMC.
#' @param log.every Integer, defining how often samples from the posterior will
#'   be written to log files and shown on screen.
# \item{monophyly}{A vector indicating monophyly constraints for the taxon sets declared with \code{taxonset}.}
# 
#  \item{tmrcaCons}{A list, containing the prior distribution(s) for age constraints of internal nodes (which must be grouped by \code{taxonset}).}
#  
#  \item{startingTree}{Either "random" or "upgma", or an object of class \code{"phylo"} to be used as a starting tree.}
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
#' # rbeauti(ips.16S, taxonset = ts)
#' @import XML
#' @export

rbeauti  <- function(..., file, template = "standard", 
                     link.clocks = TRUE, link.trees = TRUE,
                     subst.model, clock, tree,
                     taxonset, 
                     chain.length = 1e+7, log.every = 1000){
  
  ## Set file names
  ## --------------
  dir <- dirname(file)
  file <- base_name <- gsub("[.]xml$", "", basename(file))
  dir <- file.path(dir, file)
  dir <- gsub(" ", "_", paste(dir, Sys.time()))
  dir <- gsub(":", "-", dir)
  dir.create(dir)
  file <- file.path(dir, file)
  file <- paste(file, c("xml", "log", "$(tree).trees"), sep = ".")
  names(file) <- c("xml", "log", "trees")
  
  if (file.exists(file["xml"])) stop("'", file["xml"], "' already exists - please choose another name")
  
  ## Select among possible models
  ## ----------------------------
  subst.model <- match.arg(subst.model, c("JC69", "HKY", "TN93", "GTR"))
  clock <- match.arg(clock, c("Strict Clock", "Relaxed Clock Exponential", 
                              "Relaxed Clock Log Normal", "Random Local Clock"))
  tree <- match.arg(tree, c("Yule", "Calibrated Yule", "Birth Death", 
                             "Fossilized Birth Death"))
  
  ## handle partitions
  ## -----------------
  s <- list(...)
  if (unique(sapply(s, class)) == "list"){
    s <- unlist(s, recursive = FALSE)
  }
  
  ## Use tip dates
  ## -------------
  if (tree == "Fossilized Birth Death"){
    
    tip_dates <- gsub("^.*_", "", rownames(s[[1]]))
    tip_dates <- paste(paste(rownames(s[[1]]), tip_dates, sep = "="), collapse = ",")
    
  } else {
    tip_dates <- NA
  }
  
  ## assemble node(s) <data>
  ## -----------------------
  if (is.null(names(s))){
    id <- paste0("gene", 1:length(data))
    names(s) <- id
  } else {
    id <- names(s)
  }
  # id <- paste0("BEAST_", id)
  # warning("'id <- bears' active")
  id <- "bears"
  names(s) <- id
  
  
  ## Collect data (eventuall in a S4 object?)
  x <- list(
    environment = environment(),
    base.name = base_name,
    file.names = as.list(file),
    log.every = log.every,
    partition = id,
    tip.dates = tip_dates,
    subst.model = subst.model,
    clock = clock,
    link.clocks = link.clocks,
    tree = tree,
    link.trees = link.trees
  )
  
  ## Global counters to construct unique IDS
  # counter <- new.env(parent = emptyenv())
  # counter$gamma <- 0
  # counter$uniform <- 0
  # counter$realParameter <- 0
  counter <- list(realParameter = 0,
                  uniform = 0,
                  gamma = 0,
                  logNormal = 0)
  
  
  data <- assembleDataNode(s)
  
  ## assemble node <state>
  ## ---------------------
  state <- assembleStateNode(x)
  
  ## assemble node <init>
  ## ---------------------
  init <- assembleInitNode(x)
  
  ## assemble node <distribution>
  ## ----------------------------
  distribution <- assembleDistributionNode(x)
  
  ## assemble <operator> nodes
  ## -------------------------
  data("operator_list", envir = environment())
  if (x$tree != "Fossilized Birth Death"){
    operator_list <- operator_list[operator_list$tip.date != is.na(x$tip.dates), ]
  }
  ops <- operator_list$name[operator_list$clock == x$clock & operator_list$tree == x$tree]
  # if (x$subst.model == "JC69") ops <- c("KappaScaler.s", ops, "FrequenciesExchanger.s")
  if (x$subst.model == "HKY") ops <- c("KappaScaler.s", ops, "FrequenciesExchanger.s")
  if (x$subst.model == "TN93") ops <- c("kappa1Scaler.s", "kappa2Scaler.s", ops, "FrequenciesExchanger.s")
  if (x$subst.model == "GTR") ops <- c("RateACScaler.s", "RateAGScaler.s", "RateATScaler.s", 
                                        "RateCGScaler.s","RateGTScaler.s", ops, "FrequenciesExchanger.s")
  operators <- lapply(ops, operator, x = x)
  
  # operators <- assembleOperators(id, link.clocks = link.clocks, link.trees = link.trees)
  
  ## assemble node <run> 
  ## ---------------------
  run <- xmlNode("run", 
                 attrs = c(id = "mcmc",
                           spec = "MCMC",
                           chainLength = format(chain.length, scientific = FALSE)))
  run <- addChildren(run, kids = list(state))
  run <- addChildren(run, kids = init)
  run <- addChildren(run, kids = list(distribution))
  run <- addChildren(run, kids = operators)
  run <- addChildren(run, kids = assembleLoggers(x))
  
  ## assemble nodes <map> 
  ## --------------------
  m <- mm <- c("Uniform", "Exponential", "LogNormalDistributionModel", "Normal", "Beta", 
               "Gamma", "LaplaceDistribution", "Prior", "InverseGamma", "OneOnX")
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
                                      required = "",
                                      version = "2.6"))
  beast <- addChildren(beast, kids = data)
  beast <- addChildren(beast, kids = maps)
  beast <- addChildren(beast, kids = list(run))
  ## convert from class XMLNode to XMLInternalDocument
  ## -------------------------------------------------
  #   beast <- saveXML(beast, encoding = "UTF-8",
  #                    prefix = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
  #   beast <- xmlInternalTreeParse(beast, asText = TRUE)
  
  
  
  saveXML(beast, file = x$file.names$xml, encoding = "UTF-8",
          prefix = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
  
  invisible(x)
}
