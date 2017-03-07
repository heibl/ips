## This code is part of the ips package
## © C. Heibl 2015 (last update 2016-11-23)

#' @export
  
mrbayes.prset <- function(...){
  
  args <- list(...)
  args <- lapply(args, format, scientific = FALSE, trim = TRUE)
  
  arg.set <- list(
    traitiopr = NA,
    revmatpr = NA,
    aamodelpr = NA,
    aarevmatpr = NA,
    omegapr = NA,
    ny98omega1pr = NA,
    ny98omega3pr = NA,
    m3omegapr = NA,
    codoncatfreqs = NA,
    statefreqpr = NA,
    shapepr = NA,
    ratecorrpr = NA,
    pinvarpr = NA,
    covswitchpr = NA,
    symdirihyperpr = NA,
    topologypr = NA,
    brlenspr = NA,
    clockvarpr = NA,
    igrvarpr = NA
  )
  
  
  no.arg <- setdiff(names(args), names(arg.set))
  if ( length(no.arg) > 0 ){
    stop(paste('"', no.arg[1],'" is not valid argument', sep = ""))
  }
#   string.args <- names(arg.set)[arg.set != "NUMERIC"]
#   string.args <- intersect(names(args), string.args)
#   for ( i in string.args ){
#     args[[i]] <- match.arg(args[[i]], arg.set[[i]])
#   }
  args
}