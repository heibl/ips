## This code is part of the ips package
## © C. Heibl 2014 (last update 2019-06-27)

#' @export

write.phy <- function(x, file, block.width = FALSE, 
                      strict = FALSE){
	
	datatype <- ifelse(is.numeric(x[1, 1]), 
                     "continuous", "nc")
	
	ntax <- nrow(x)
	nchar <- ncol(x)
  
	# taxon names
	# -----------
	taxnames <- rownames(x)
	
	if (strict){
		taxnames <- substring(taxnames, 1, truncate)
		missing <- 10 - unlist(lapply(strsplit(taxnames, ""), length))
		for (i in seq(along = taxnames))
			taxnames[i] <- paste(taxnames[i], paste(rep("*",
			    missing[i]), collapse = ""), sep = "")
		if (any(duplicated(taxnames))) 
		    warning("truncation of taxon names created identical strings")
	}
	else {
		xx <- nchar(taxnames)
		diff <- max(xx) - xx + 3 # add 3 for baseml
		for (i in 1:ntax) taxnames[i] <- paste(taxnames[i], 
			paste(rep(" ", diff[i]), collapse = ""), sep = "")
	}
	
	# indices of partitions: pt
	# -------------------------
	if (!block.width){
    block.width <- nchar
	} else {
    if ( block.width > nchar ){
      block.width <- nchar
    }
	}
	nbpart <- ceiling(nchar/block.width)
	pt <- matrix(nrow = nbpart, ncol = 2)
	pt[1, ] <- c(1, block.width)
	if (nbpart > 1)
		for (i in 2:(dim(pt)[1])){
		    pt[i, ] <- c(pt[i - 1, 2] + 1, pt[i - 1, 2] + block.width)
		    pt[nbpart, 2] <- nchar
		}
	
	# assemble matrix: m
	# ------------------
	phy <- paste(ntax, nchar)
	for (i in seq(along = pt[, 1])){
		sm <- as.character(x[, pt[i, 1]:pt[i, 2]])
		if (is.null(dim(sm))) sm <- as.matrix(sm, ncol = 1)
		sm <- apply(sm, 1, paste, collapse = "")
		if (i == 1) sm <- paste(taxnames, sm)
		if (i < max(seq(along = pt[, 1]))) sm <- c(sm, "")
		phy <- c(phy, sm)
	}
	
	# write PHYLIP file
	# -----------------
	if (missing(file)) {
	  ## return character vector:
	  return(phy) 
	} else {
	  if ( file == "" ) {
	    ## print onto screen:
	    cat(phy, sep = "\n")
	  } else {
	    ## write to file:
	    write(phy, file = file)
	  }
	}
}
