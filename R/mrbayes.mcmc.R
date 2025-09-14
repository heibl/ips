## This code is part of the ips package
## © C. Heibl 2015 (last update 2016-11-23)

#' @export

#' @title MrBayes MCMC settings
#' @description
#' Specify the MCMC settings for \code{\link{mrbayes}}.
#' @param ngen A character string
#' @param nruns xxx
#' @param nchains xxx
#' @param temp xxx
#' @param reweight xxx
#' @param swapfreq xxx
#' @param nswaps xxx
#' @param samplefreq xxx
#' @param printfreq xxx
#' @param printall xxx
#' @param printmax xxx
#' @param mcmcdiagn xxx
#' @param diagnfreq xx
#' @param diagnstat xxx
#' @param minpartfreq xxx
#' @param allchains xxx
#' @param allcomps xxx
#' @param relburnin xxx
#' @param burnin xxx
#' @param burninfrac xxx
#' @param stoprule xxx
#' @param stopval xxx
#' @param savetrees xxx
#' @param checkpoint xxx
#' @param checkfreq xxx
#' @param startparams xxx
#' @param starttree xxx
#' @param nperts xxx
#' @param data xxx
#' @param ordertaxa xxx
#' @param append xxx
#' @param autotune xxx
#' @param tunefreq xxx
#' @details
#' xxx
#' @returns A list of MCMC parameters.
#' @seealso \code{\link{mrbayes}}, \code{\link{mrbayes.lset}} , \code{\link{mrbayes.prset}}.
#' @export

mrbayes.mcmc <- function(ngen = 1000000, nruns = 2, nchains = 4, temp = 0.1,
                         reweight = c(0, 0), swapfreq = 1, nswaps = 1, 
                         samplefreq = 500, printfreq = 1000, printall = TRUE,
                         printmax = 8, mcmcdiagn = TRUE, diagnfreq = 5000, 
                         diagnstat = "avgstddev", minpartfreq = 0.1, 
                         allchains = FALSE, allcomps = FALSE, relburnin = TRUE,
                         burnin = 0, burninfrac = 0.25, stoprule = FALSE, 
                         stopval = 0.05, savetrees = FALSE, checkpoint = TRUE,
                         checkfreq = 2000, startparams = "current", 
                         starttree = "current", nperts = 0, data = TRUE,
                         ordertaxa = FALSE, append = FALSE, autotune = TRUE,
                         tunefreq = 100){
  
  arg.set <- list(
    ngen = "NUMERIC",
    nruns = "NUMERIC",
    nchains = "NUMERIC",
    temp = "NUMERIC",
    # reweight        <number>,<number>     0.00 v 0.00 ^,
    swapfreq = "NUMERIC",                                      
    nswaps = "NUMERIC",                                     
    samplefreq = "NUMERIC",                                      
    printfreq = "NUMERIC",
    printall = c("yes", "no"),
    printmax = "NUMERIC",
    mcmcdiagn = c("yes", "no"),
    diagnfreq = "NUMERIC",
    diagnstat = c("avgstddev", "maxstddev"),
    minpartfreq = "NUMERIC",
    allchains = c("yes", "no"),
    allcomps = c("yes", "no"),
    relburnin = c("yes", "no"),
    burnin = "NUMERIC",                                     
    burninfrac = "NUMERIC",   
    stoprule = c("yes", "no"),
    stopval = "NUMERIC",
    savetrees = c("yes", "no"),
    checkpoint = c("yes", "no"),
    checkfreq = "NUMERIC",
    # filename = filename,
    startparams = c("current", "reset"),
    starttree = c("current", "random", "parsimony"),
    nperts = "NUMERIC", 
    data = c("yes", "no"),
    ordertaxa = c("yes", "no"),
    append = c("yes", "no"),
    autotune = c("yes", "no"),
    tunefreq = "NUMERIC" 
  )
  
  args <- list(ngen = ngen, nruns = nruns, nchains = nchains, temp = temp, 
               reweight = reweight, swapfreq = samplefreq, nswaps = nswaps, 
               samplefreq = samplefreq, printfreq = printfreq, 
               printall = printall, printmax = printmax, mcmcdiagn = mcmcdiagn, 
               diagnfreq = diagnfreq, diagnstat = diagnstat, 
               minpartfreq = minpartfreq, allchains = allchains, 
               allcomps = allcomps, relburnin = relburnin, burnin = burnin, 
               burninfrac = burninfrac, stoprule = stoprule, stopval = stopval, 
               savetrees = savetrees, checkpoint = checkpoint, 
               checkfreq = checkfreq, startparams = startparams, 
               starttree = starttree, nperts = nperts, data = data, 
               ordertaxa = ordertaxa, append = append, autotune = autotune,
               tunefreq = tunefreq)
  
  # no.arg <- setdiff(names(args), names(arg.set))
  # if ( length(no.arg) > 0 ){
  #   stop(paste('"', no.arg[1],'" is not valid argument', sep = ""))
  # }
  # string.args <- names(arg.set)[arg.set != "NUMERIC"]
  # string.args <- intersect(names(args), string.args)
  # for ( i in string.args ){
  #   args[[i]] <- match.arg(args[[i]], arg.set[[i]])
  # }
  args
}