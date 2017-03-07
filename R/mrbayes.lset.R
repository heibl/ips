## This code is part of the ips package
## © C. Heibl 2015 (last update 2015-11-28)

#' @export

mrbayes.lset <- function(..., partition){
  
  args <- list(...)
  args <- lapply(args, format, scientific = FALSE, trim = TRUE)
  
  arg.set <- list(
    nucmodel = c("4by4", "doublet", "codon", "protein"), 
    nst = c("1", "2", "6", "mixed"), 
    code = c("universal", 
             "vertmt", 
             "mycoplasma", 
             "yeast", 
             "ciliates", 
             "metmt"), 
    ploidy = c("haploid", "diploid", "zlinked"), 
    rates = c("equal", 
              "gamma", 
              "propinv", 
              "invgamma",
              "adgamma"), 
    ngammacat = format(1:24, trim = TRUE), 
    nbetacat = format(1:24, trim = TRUE),  
    omegavar = c("equal", "ny98", "m3"),
    covarion = c("no", "yes"), 
    coding = c("all", 
               "variable", 
               "noabsencesites",
               "nopresencesites"), 
    parsmodel = c("no", "yes")
  )
  no.arg <- setdiff(names(args), names(arg.set))
  if ( length(no.arg) > 0 ){
    stop(paste('"', no.arg[1],'" is not valid argument', sep = ""))
  }
  for ( i in names(args) ){
    args[[i]] <- match.arg(args[[i]], arg.set[[i]])
  }
  
  ## partition
  if ( !missing(partition) ){
    attr(args, "partition") <- partition
  }
  
  args
}