## This code is part of the ips package
## © C. Heibl 2015 (last update 2025-09-13)

#' @title MrBayes model settings
#' @description
#' Specify a model of DNA, Protein or trait evolution for \code{\link{mrbayes}}.
#' @param nucmodel A character string
#' @param nst xxx
#' @param code xxx
#' @param ploidy xxx
#' @param rates xxx
#' @param ngammacat xxx
#' @param nlnormcat xxx
#' @param nmixtcat xxx
#' @param nbetacat xxx
#' @param omegavar xxx
#' @param covarion xxx
#' @param coding xxx
#' @param parsmodel xxx
#' @param partition A character string serving as an identifier in a
#'   multi-partition model.
#' @details
#' xxx
#' @returns A list of model parameters.
#' @seealso \code{\link{mrbayes}}, \code{\link{mrbayes.prset}} , \code{\link{mrbayes.mcmc}}.
#' @export

mrbayes.lset <- function(nucmodel = "4by4", nst = 1, code = "universal", 
                         ploidy = "diploid", rates = "equal", ngammacat = 4,
                         nlnormcat = 4, nmixtcat = 4, nbetacat = 5, 
                         omegavar = "equal", covarion = "no", coding = "all",
                         parsmodel = "no", partition = 1){
  
  
  arg_set <- list(
    nucmodel = c("4by4", "doublet", "codon", "protein"),
    nst = c("1", "2", "6", "mixed"),
    code = c("universal", "vertmt", "mycoplasma", "yeast", "ciliates", "metmt"),
    ploidy = c("haploid", "diploid", "zlinked"),
    rates = c("equal", "gamma", "propinv", "invgamma", "adgamma"),
    ngammacat = format(1:24, trim = TRUE),
    nlnormcat = format(1:24, trim = TRUE),
    nmixtcat = format(1:24, trim = TRUE),
    nbetacat = format(1:24, trim = TRUE),
    omegavar = c("equal", "ny98", "m3"),
    covarion = c("no", "yes"),
    coding = c("all", "variable", "noabsencesites", "nopresencesites"),
    parsmodel = c("no", "yes")
  )
  
  ## Dieser Code ist zu flexibel für den vorliegenden Fall
  # args <- list(...)
  # args <- lapply(args, format, scientific = FALSE, trim = TRUE)
  # 
  # arg.set <- list(
  #   nucmodel = c("4by4", "doublet", "codon", "protein"), 
  #   nst = c("1", "2", "6", "mixed"), 
  #   code = c("universal", 
  #            "vertmt", 
  #            "mycoplasma", 
  #            "yeast", 
  #            "ciliates", 
  #            "metmt"), 
  #   ploidy = c("haploid", "diploid", "zlinked"), 
  #   rates = c("equal", 
  #             "gamma", 
  #             "propinv", 
  #             "invgamma",
  #             "adgamma"), 
  #   ngammacat = format(1:24, trim = TRUE), 
  #   nbetacat = format(1:24, trim = TRUE),  
  #   omegavar = c("equal", "ny98", "m3"),
  #   covarion = c("no", "yes"), 
  #   coding = c("all", 
  #              "variable", 
  #              "noabsencesites",
  #              "nopresencesites"), 
  #   parsmodel = c("no", "yes")
  # )
  # args <- args[match(names(args), names(arg.set))]
  # no.arg <- setdiff(names(args), names(arg.set))
  # if (length(no.arg) > 0){
  #   stop(paste('"', no.arg[1],'" is not valid argument', sep = ""))
  # }
  # for (i in names(args)){
  #   args[[i]] <- match.arg(args[[i]], arg.set[[i]])
  # }
  
  args <- list(nucmodel, nst, code, ploidy, rates, ngammacat, nlnormcat, 
               nmixtcat, nbetacat, omegavar, covarion, coding, parsmodel)
  names(args) <- names(arg_set)
  ## partition
  if (!missing(partition)){
    attr(args, "partition") <- partition
  }
  
  args
}