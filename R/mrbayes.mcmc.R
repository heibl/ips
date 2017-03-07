## This code is part of the ips package
## © C. Heibl 2015 (last update 2016-11-23)

#' @export

mrbayes.mcmc <- function(...){
  
  args <- list(...)
  args <- lapply(args, format, scientific = FALSE, trim = TRUE)
  
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
  
  
  no.arg <- setdiff(names(args), names(arg.set))
  if ( length(no.arg) > 0 ){
    stop(paste('"', no.arg[1],'" is not valid argument', sep = ""))
  }
  string.args <- names(arg.set)[arg.set != "NUMERIC"]
  string.args <- intersect(names(args), string.args)
  for ( i in string.args ){
    args[[i]] <- match.arg(args[[i]], arg.set[[i]])
  }
  args
}