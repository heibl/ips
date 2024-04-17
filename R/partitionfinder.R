## This code is part of the megaptera package
## © C. Heibl 2014 (last update 2016-11-08)

#' @title PartitionFinder
#' @description Provides a wrapper to the PartitionFinder software.
#' @param alignment A
#' @param user.tree A
#' @param branchlengths A
#' @param models A
#' @param model.selection A
#' @param search A
#' @param exec A character string giving the path to the executable (python 
#'   script).
#' @references 
#'   Lanfear, R., B. Calcott, S.Y.W. Ho, and S. Guindon. 2012. 
#'   PartitionFinder: combined selection of partitioning schemes and
#'   substitution models for phylogenetic analyses. \emph{Molecular 
#'   Biology and Evolution} \strong{29}: 1695-1701.  
#'   
#'   Lanfear, R., B. Calcott, K. David, C. Mayer, and A. Stamatakis. 2014. 
#'   Selecting optimal partitioning schemes for phylogenomic datasets.
#'   \emph{BMC Evolutionary Biology} \strong{14}: 82.
#' @export

partitionfinder <- function(alignment, user.tree, 
                            branchlengths = "linked", models = "all",
                            model.selection = "BIC", search = "greedy",
  exec = "/Applications/PartitionFinderV1.1.1_Mac/PartitionFinder.py"){
  
  ## check arguments
  branchlengths <- match.arg(branchlengths, c("linked", "unlinked"))
  ## extension pending: <list>
  models <- match.arg(models, c("all", "raxml", "mrbayes", "beast"))
  model.selection <- match.arg(models, c("AIC", "AICc", "BIC"))
  
  write.phy(alignment, file = "test.phy")
  
  ## data blocks
  ## -----------
  # data.blocks
  
  ## create configuration file
  ## -------------------------
  cfg <- c(
    "alignment = test.phy",
    paste0("branchlengths = ", branchlengths, ";"),
    paste0("models = ", models, ";"),
    paste0("model_selection = ", model.selection, ";"),
    "[data_blocks]",
    paste0("search = ", search, ";")
  )
  write(cfg, file = "partition_finder.cfg")
  
  
  system(paste("python",
               exec,
               "/Applications/PartitionFinderV1.1.1_Mac/examples/nucleotide"))
}
