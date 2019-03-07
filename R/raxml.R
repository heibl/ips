## This code is part of the ips package
## © C. Heibl 2014 (last update 2019-02-05)

#'@title Maximum Likelihood Tree Estimation with RAxML
#'@description Provides an interface to the C program \bold{RAxML} (see
#'  Reference section) for maximum likelihood estimation of tree topology and/or
#'  branch lengths, rapid and conventional non-parametric bootstrapping, mapping
#'  splits onto individual topologies, and a lot more. See the RAxML manual for
#'  details, especially if you are a new user of RAxML.
#'@param DNAbin A matrix of DNA sequences of class \code{\link{DNAbin}}.
#'@param m A vector of mode \code{"character"} defining a model of molecular
#'  evolution; currently only GTR model available.
#'@param f A vector of mode \code{"character"} selecting an RAxML algorithm
#'  analogous to the \code{-f} flag (see Detail section and RAxML manual).
#'@param N Either of mode \code{"integer"} or \code{"character"}. Integers give
#'  the number of independant searches on different starting tree or replicates
#'  in bootstrapping. Alternatively, one of four bootstopping criteria can be
#'  chosen: \code{"autoFC"}, \code{"autoMR"}, \code{"autoMRE"}, or
#'  \code{"autoMRE_IGN"}.
#'@param p Integer, setting a random seed for the parsimony starting trees.
#'@param b Integer, setting a random seed for bootstrapping.
#'@param x Integer, setting a random seed for rapid bootstrapping.
#'@param k Logical, if \code{TRUE}, the branch lengths of bootstrapped trees are
#'  recorded.
#'@param weights A vector of mode \code{"numeric"} giving integers to assign
#'  individual weights to each column of the alignment. (-a)
#'@param partitions A data frame giving the partitions of the alignment.
#'@param outgroup A vector of mode \code{"character"} containing the names of
#'  the outgroup taxa.
#'@param backbone A \code{\link{phylo}} object representing a backbone tree.
#'@param file A vector of mode \code{"character"} giving a name for the output
#'  files.
#'@param exec A vector of mode \code{"character"} giving the path to the
#'  directory containing the RAxML executable. The default value will work on
#'  Mac OS X if the folder containing the executale is renamed to
#'  \code{"RAxML-8.0.3"}.
#'@param threads Integer, giving the number of parallel threads to use (PTHREADS
#'  only).
#'@details There are some limitations of this wrapper compared to RAxML run
#'  directly from the command line. \enumerate{ \item Only DNA is allowed as
#'  data type. \item Option \code{f} can only take a limited number of values
#'  (\code{d}, \code{a}). } % close enumerate
#'
#'  RAxML needs the specification of random seeds for parsimony estimation of
#'  starting trees and for bootstrap resampling. The corresponding argument
#'  names in \code{raxml} are identical to the flags used by RAxML (\code{-p},
#'  \code{-b}, and \code{-x}). If you choose not to give any values,
#'  \code{raxml} will generate a (different) value for each requiered random
#'  seed every time it is called. Be aware that \code{\link{set.seed}} will work
#'  only for \code{p}, but not for \code{b} or \code{x}.
#'@return A list with a variable number of elements, depending on the analysis
#'  chosen: \tabular{ll}{ \code{"info"} \tab RAxML log file as character
#'  string\cr \code{"bestTree"} \tab MLE of tree\cr \code{"bipartitions"} \tab
#'  MLE of tree annotated with bootstrap proportions\cr \code{"bootstrap"} \tab
#'  bootstrapped trees\cr }
#'@references (in chronolocigal order)
#'
#'  Stamatakis, A., T. Ludwig and H. Meier. 2004. RAxML-III: A fast program for
#'  maximum likelihood-based inference of large phylogenetic trees.
#'  \emph{Bioinformatics} \bold{1}: 1--8.
#'
#'  Stamatakis, A. 2006. RAxML-VI-HPC: Maximum likelihood-based phylogenetic
#'  analyses with thousands of taxa and mixed models. \emph{Bioinformatics}
#'  \bold{22}: 2688--2690.
#'
#'  Stamatakis, A., P. Hoover, and J. Rougemont. 2008. A rapid bootstrap
#'  algorithm for the RAxML web-servers. \emph{Syst. Biol.} \bold{75}: 758--771.
#'
#'  Pattengale, N. D., M. Alipour, O. R. P. Bininda-Emonds, B. M. E. Moret, and
#'  A. Stamatakis. 2010. How many bootstrap replicates are necessary?
#'  \emph{Journal of Computational Biology} \bold{17}: 337-354.
#'
#'  Stamatakis, A. 2014. RAxML Version 8: A tool for phylogenetic analysis and
#'  post-analysis of large phylogenies. \emph{Bioinformatics} Advance Access.
#'@note RAxML is a C program and the source code is not contained in this
#'  package. This means that in order to run this function you will need to
#'  install RAxML yourself. See the 'Software' tab on
#'  \url{http://www.exelixis-lab.org/} for the most recent documentation and
#'  source code of RAxML. Depending on where you chose to install RAxML, you
#'  need to adjust the \code{exec} argument.
#'
#'  \code{raxml} was last tested and running fine on Mac OS X with RAxML 8.0.29.
#'  Please be aware that calling third-party software from within R is a
#'  platform-specific process and I cannot guarantee that \code{raxml} will
#'  behave properly on any system.
#'@seealso \code{\link{raxml.partitions}} to store partitioning information in a
#'  data frame suitable for input as \code{partitions} argument in \code{raxml}.
#'@examples
#'## bark beetle sequences
#'data(ips.cox1)
#'data(ips.16S)
#'data(ips.28S)
#'
#'ips <- cbind(ips.cox1, ips.16S, ips.28S,
#'             fill.with.gaps = TRUE)
#'
#'exec <- "/Applications/RAxML-code/standard-RAxML/raxmlHPC-PTHREADS-AVX"
#'w <- sample(1:5, ncol(ips.cox1), replace = TRUE)
#'
#'\dontrun{
#'
#'# Simple tree search with GTRCAT and GTRGAMMA
#'tr <- raxml(ips.cox1, f = "d", N = 2, p = 1234,
#'            exec = exec) # -1743.528461
#'tr <- raxml(ips.cox1, m = "GTRGAMMA", f = "d", N = 2, p = 1234,
#'            exec = exec)
#'    
#'# Applying weights to columns                   
#'tr <- raxml(ips.cox1, f = "d", N = 2, p = 1234,
#'            weights = w, exec = exec) # -1743.528461
#'
#'# Rapid bootstrap
#'tr <- raxml(ips.cox1, m = "GTRGAMMA",
#'            f = "a", N = 10, p = 1234, x = 1234,
#'            exec = exec)
#'
#'# Rapid bootstrap with automatic halt
#'tr <- raxml(ips.cox1, m = "GTRGAMMA",
#'            f = "a", N = "autoMRE", p = 1234, x = 1234,
#'            exec = exec)
#'}
#'@export

raxml <- function(DNAbin, m = "GTRCAT", f, N, p, b, x, k,
                  weights, partitions, outgroup, backbone = NULL, 
                  file = "fromR", exec, threads){
  
  if (!inherits(DNAbin, "DNAbin")) stop("DNAbin is not of class 'DNAbin'")
				
	# number of threads (PTHREADS only)
	# ---------------------------------
	if (!missing(threads))
    exec <- paste(exec, "-T", threads)
	
	## input file names
	## ----------------
	rin <- c(s = paste("-s ", file, ".phy", sep = ""),
	         n = paste("-n", file),
	         napt = paste("-n ", file, ".APT", sep = ""),
	         weight = paste0(file, "-weights.txt"),
	         partition = paste0(file, "-partitions.txt"),
	         tree = paste0(file, "-backbone.tre"))
  
	# Clear previous runs with same tag
	# ---------------------------------
	unlink(rin)
	unlink(list.files(pattern = paste0("RAxML_[[:alpha:]]+[.]", file))) 
  
  # Substitution model
  # ------------------
	m <- match.arg(m, c("GTRCAT", "GTRCATX", 
	                    "GTRCATI", "GTRCATIX",
	                    "ASC_GTRCAT", "ASC_GTRCATX", 
	                    "GTRGAMMA", "GTRGAMMAX", 
	                    "GTRGAMMAI", "GTRGAMMAIX",
	                    "ASC_GTRGAMMA", "ASC_GTRGAMMAX"))
  m <- paste("-m", m)
  
  ## Number of searches/replicates
  ## -----------------------------
  if (is.character(N)){
    N <- match.arg(N, c("autoFC", "autoMR", "autoMRE", "autoMRE_IGN"))
  }
  N <- paste("-N", N)
  
  ## Random seeds
  ## ------------
	rs <- function(rseed, type = "p"){
	  if (missing(rseed)) rseed <- sample(1:999999, 1)
	  paste0("-", type, " ",  rseed)
	}
	p <- rs(p)
	if (!missing(x)) x <- rs(x, type = "x")
  
	## Write sequences to input file
	## -----------------------------
	write.phy(DNAbin, paste(file, "phy", sep = "."))

  ## rout: raxml output file names
  ## -----------------------------
  output.types <- c("info", "bestTree", "bootstrap", "bipartitions")
  rout <- paste("RAxML_", output.types, ".", file, sep = "")
  names(rout) <- output.types
  
  ## Algorithms
  ## ----------
  if (missing(f)) f <- "d"
  f <- match.arg(f, c("d", "a"))
  alg <- paste("-f", f, p) # add parsimony seed to algorithm
  if (f == "a") alg <- paste(alg, x)
  if (!missing(b)) alg <- paste(alg, "-b", b)
  if (missing(N)) stop("the number of runs must be given (N)")
  
  ## Outgroup
  ## --------
  if (missing(outgroup)){
    o <- ""
  } else {
    if (length(grep(",", outgroup))) outgroup <- unlist(strsplit(outgroup, ","))
    o <- outgroup %in% rownames(DNAbin)
    if (!all(o)){
      o <- paste(paste("\n  -", outgroup[!o]), collapse = "")
      stop(paste("outgroup names not in 'DNAbin':",   o))
    }
    o <- paste(outgroup, collapse = ",")
    o <- paste("-o", o)
  }
  
  ## Columns weights
  ## ---------------
  if (!missing(weights)){
    if (length(weights) != ncol(DNAbin)) stop("number of 'weights' don't equal number of columns")
    write(weights, rin["weight"], length(weights))
    weights <- paste("-a", rin["weight"])
  } else {
    weights <- ""
  }

	# Write partition file
	## ------------------
	if (!missing(partitions)) {
    if (is.character(partitions)){
      q <- partitions
    } else {
      q <- paste(partitions$type, ", ", 
                 partitions$locus, " = ", 
                 partitions$begin, "-", 
                 partitions$end, sep = "")
    }
		write(q, rin["partition"])
    multipleModelFileName <- paste(" -q", rin["partition"], "")
	} else multipleModelFileName <- ""

  ## Backbone tree
  ## -------------
	if (!is.null(backbone)){
	  write.tree(backbone, rin["tree"])
	  g <- paste(" -g", rin["tree"], "")
	} else {
	  g <- " "
	}
  
  ## Save branch lengths of bootstrap replicates
  ## -------------------------------------------
  if (missing(k)) k <- FALSE 
  k <- ifelse(k, "-k", "")
  
	## Prepare and execute call
  ## ------------------------
	CALL <- paste(exec, alg, m, o, k,
	              weights, 
	              multipleModelFileName, N, g, 
	              rin["s"], rin["n"])
  
	if (length(grep("MPI", exec))) system(paste("mpirun", CALL))
  print(CALL)
	system(CALL)
  
	res <- scan(rout["info"], quiet = TRUE, what = "char", sep = "\n")
	if (length(grep("exiting", res)))
	  stop("\n", paste(res, collapse = "\n"))
	
	## Read results
	## ------------
	bestTree <- bipartitions <- bootstrap <- NULL
  info <- scan(rout["info"], what = "c", sep = "\n", quiet = TRUE)
	if (f %in% c("a", "d") & missing(b))
	  bestTree <- read.tree(rout["bestTree"])
	if (f %in% c("a"))
	  bipartitions <- read.tree(rout["bipartitions"])
	if (!missing(b) | !missing(x))
	  bootstrap <- read.tree(rout["bootstrap"])
	
	obj <- list(info = info,
	            bestTree = bestTree,
	            bipartitions = bipartitions,
	            bootstrap = bootstrap)
	obj[sapply(obj, is.null)] <- NULL
	obj
}

# ##########################################################
# # 	test initial rearrangement
# ##########################################################
# 
# if ( optimize ){
#   n <- 0:4
#   call.opt <- c(paste("./raxmlHPC ", rs(), " -y -s ", file , " -m GTRCAT -n ST", n, sep = ""),
#                 paste("./raxmlHPC -f d -i 10 -m GTRCAT -s ", file ,
#                       " -t RAxML_parsimonyTree.ST", n, " -n FI", n,  sep = ""),
#                 paste("./raxmlHPC -f d -m GTRCAT -s ", file , 
#                       " -t RAxML_parsimonyTree.ST", n, " -n AI", n, sep = ""))
#   lapply(call.opt, system)
#   fixed <- paste("FI", n, sep = "")
#   auto <- paste("AI", n, sep = "")
#   res <- paste("RAxML_info", c(fixed, auto), sep = ".")
#   res <- lapply(res, scan, what = "c", sep = "", quiet = TRUE)
#   res <- as.numeric(sapply(res, function(x) x[grep("Program", x) - 1]))
#   names(res) <- c(fixed, auto)
#   FIXED <- mean(res[fixed])
#   AUTO <- mean(res[auto])
#   # set i
#   if ( AUTO > FIXED ) {
#     i <- x[grep("setting", x) + 1]
#     i <- as.numeric(gsub(",", "", i))
#   }
#   
#   
#   ##########################################################
#   # 	number of categories
#   ##########################################################
#   x <- paste("./raxmlHPC -f d -i ", i, " -m GTRCAT -s ", file ," -t RAxML_parsimonyTree.ST", sep = "")
#   system(paste(x, "0 -c 10 -n C10_0", sep = ""))
#   system(paste(x, "1 -c 10 -n C10_1", sep = ""))
#   system(paste(x, "2 -c 10 -n C10_2", sep = ""))
#   system(paste(x, "3 -c 10 -n C10_3", sep = ""))
#   system(paste(x, "4 -c 10 -n C10_4", sep = ""))
#   system(paste(x, "0 -c 40 -n C40_0", sep = ""))
#   system(paste(x, "1 -c 40 -n C40_1", sep = ""))
#   system(paste(x, "2 -c 40 -n C40_2", sep = ""))
#   system(paste(x, "3 -c 40 -n C40_3", sep = ""))
#   system(paste(x, "4 -c 40 -n C40_4", sep = ""))
#   system(paste(x, "0 -c 55 -n C55_0", sep = ""))
#   system(paste(x, "1 -c 55 -n C55_1", sep = ""))
#   system(paste(x, "2 -c 55 -n C55_2", sep = ""))
#   system(paste(x, "3 -c 55 -n C55_3", sep = ""))
#   system(paste(x, "4 -c 55 -n C55_4", sep = ""))
#   x <- scan("RAxML_info.C10_0", what="c", sep = "", 
#             quiet = TRUE)
#   C10_0 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C10_1", what="c", sep = "", 
#             quiet = TRUE)
#   C10_1 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C10_2", what="c", sep = "", 
#             quiet = TRUE)
#   C10_2 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C10_3", what="c", sep = "", 
#             quiet = TRUE)
#   C10_3 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C10_4", what="c", sep = "", 
#             quiet = TRUE)
#   C10_4 <- x[grep("Final", x)-1]
#   C10 <- mean(as.numeric(c(C10_0, C10_1, C10_2, C10_3, 			C10_4)))
#   
#   x <- scan("RAxML_info.C40_0", what="c", sep = "", quiet = TRUE)
#   C40_0 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C40_1", what="c", sep = "", quiet = TRUE)
#   C40_1 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C40_2", what="c", sep = "", quiet = TRUE)
#   C40_2 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C40_3", what="c", sep = "", quiet = TRUE)
#   C40_3 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C40_4", what="c", sep = "", quiet = TRUE)
#   C40_4 <- x[grep("Final", x)-1]
#   C40 <- mean(as.numeric(c(C40_0, C40_1, C40_2, C40_3, C40_4)))
#   
#   x <- scan("RAxML_info.C55_0", what = "c", sep = "", 		quiet = TRUE)
#   C55_0 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C55_1", what = "c", sep = "", 		quiet = TRUE)
#   C55_1 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C55_2", what = "c", sep = "", 		quiet = TRUE)
#   C55_2 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C55_3", what = "c", sep = "", 		quiet = TRUE)
#   C55_3 <- x[grep("Final", x)-1]
#   x <- scan("RAxML_info.C55_4", what = "c", sep = "", 		quiet = TRUE)
#   C55_4 <- x[grep("Final", x)-1]
#   C55 <- mean(as.numeric(c(C55_0, C55_1, C55_2, C55_3, C55_4)))
#   
#   C <- if (FIXED > AUTO) c(C10, FIXED, C40, C55) 		else c(C10, AUTO, C40, C55)
#   catnum <- c(10, 25, 40, 55)
#   DF <- cbind(catnum, C)
#   DF <- DF[order(DF[,2], decreasing = TRUE),]
#   numcat <- DF[1,1]
#   
# } # end of OPTIMIZE
# else {
#   i <- 10
#   numcat = 25
# }
