## This code is part of the ips package
## Written by C. Heibl 2015 (last update 2025-09-14)

#' @title MrBayes prior settings
#' @description Specify parameters for prior distributions for
#' \code{\link{mrbayes}}.
#' @param traitiopr A character string
#' @param revmatpr xxx
#' @param aamodelpr xxx
#' @param omegapr xxx
#' @param aarevmatpr xxx
#' @param omegapr xxx
#' @param ny98omega1pr xxx
#' @param ny98omega3pr xxx
#' @param m3omegapr xxx
#' @param codoncatfreqs xxx
#' @param statefreqpr xxx
#' @param shapepr xxx
#' @param ratecorrpr xxx
#' @param pinvarpr xxx
#' @param covswitchpr xxx
#' @param symdirihyperpr xxx
#' @param topologypr xxx
#' @param brlenspr xxx
#' @param treeagepr xxx
#' @param speciationpr xxx
#' @param extinctionpr xxx
#' @param fossilizationpr xxx
#' @param samplestrat xxx
#' @param sampleprop xxx
#' @param popsizepr xxx
#' @param popvarpr xxx
#' @param nodeagepr xxx
#' @param clockratepr xxx
#' @param clockvarpr xxx
#' @param cppratepr xxx
#' @param cppmultdevpr xxx
#' @param tk02varpr xxx
#' @param wnvarpr xxx
#' @param igrvarpr xxx
#' @param ilnvarpr xxx
#' @param ratepr xxx
#' @param generatepr xxx
#' @details xxx
#' @returns A list of parameters for prior distributions.
#' @seealso \code{\link{mrbayes}}, \code{\link{mrbayes.lset}} ,
#'   \code{\link{mrbayes.mcmc}}.
#' @export

mrbayes.prset <- function(traitiopr, revmatpr, aamodelpr, aarevmatpr, omegapr, 
                          ny98omega1pr, 
                          ny98omega3pr, m3omegapr, codoncatfreqs, statefreqpr, 
                          shapepr, ratecorrpr, pinvarpr, covswitchpr, 
                          symdirihyperpr, topologypr, brlenspr, treeagepr, 
                          speciationpr, extinctionpr, fossilizationpr, 
                          samplestrat, sampleprop, popsizepr, popvarpr,
                          nodeagepr, clockratepr, clockvarpr, cppratepr, 
                          cppmultdevpr, tk02varpr, wnvarpr, igrvarpr, ilnvarpr, 
                          ratepr, generatepr){

  # args <- lapply(args, format, scientific = FALSE, trim = TRUE)
  
  # arg_set <- list(
  #   traitiopr, revmatpr, aamodelpr, omegapr, ny98omega1pr, ny98omega3pr,
  #   m3omegapr, codoncatfreqs, statefreqpr, shapepr, ratecorrpr, pinvarpr,
  #   covswitchpr, symdirihyperpr, topologypr, brlenspr, treeagepr, speciationpr,
  #   extinctionpr, fossilizationpr, samplestrat, sampleprop, popsizepr, popvarpr,
  #   nodeagepr, clockratepr, clockvarpr, cppratepr, cppmultdevpr, tk02varpr, 
  #   wnvarpr, igrvarpr, ILNvarpr, ratepr, generatepr
  # )
  
  args <- list(
    traitiopr = traitiopr, revmatpr = traitiopr, aamodelpr = aamodelpr, 
    omegapr = omegapr, ny98omega1pr = ny98omega1pr, ny98omega3pr = ny98omega3pr,
    m3omegapr = m3omegapr, codoncatfreqs = codoncatfreqs, 
    statefreqpr = statefreqpr, shapepr = shapepr, ratecorrpr = ratecorrpr, 
    pinvarpr = pinvarpr, covswitchpr = covswitchpr, 
    symdirihyperpr = symdirihyperpr, topologypr = topologypr, brlenspr = brlenspr, 
    treeagepr = treeagepr, speciationpr = speciationpr, 
    extinctionpr = extinctionpr, fossilizationpr = fossilizationpr, 
    samplestrat = samplestrat, sampleprop = sampleprop, popsizepr = popsizepr, 
    popvarpr = popvarpr, nodeagepr = nodeagepr, clockratepr = clockratepr, 
    clockvarpr = clockvarpr, cppratepr = cppratepr, cppmultdevpr = cppmultdevpr, 
    tk02varpr = tk02varpr, wnvarpr = wnvarpr, igrvarpr = igrvarpr, 
    ilnvarpr = ilnvarpr, ratepr = ratepr, generatepr = generatepr
  )
  
  
  # no.arg <- setdiff(names(args), names(arg.set))
  # if ( length(no.arg) > 0 ){
  #   stop(paste('"', no.arg[1],'" is not valid argument', sep = ""))
  # }
#   string.args <- names(arg.set)[arg.set != "NUMERIC"]
#   string.args <- intersect(names(args), string.args)
#   for ( i in string.args ){
#     args[[i]] <- match.arg(args[[i]], arg.set[[i]])
#   }
  args
}