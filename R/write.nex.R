## This code is part of the ips package
## © C. Heibl 2014 (last update 2015-12-10)

write.nex <- function(x, file, block.width = 60, 
                      taxblock = FALSE){
  
  if ( !is.list(x) ) x <- list(x)
  
  ## data types of partitions
  datatype <- sapply(x, class)
  datatype[datatype == "DNAbin"] <- "dna"
  datatype[datatype == "dist"] <- "distances"
  datatype[datatype == "data.frame"] <- "standard"
  m <- function(x) ifelse(as.raw(2) %in% x, "?", "N")
  p <- cbind(rep(1, length(x)), sapply(x, ncol))
  if ( nrow(p) > 1 ) {
    for ( i in 2:nrow(p) ){
      p[i, 1] <- p[i - 1, 2] + 1
      p[i, 2] <- p[i, 1] + p[i, 2] -1
    }
  }
  info <- data.frame(datatype = datatype,
                     missing = sapply(x, m),
                     length = sapply(x, ncol),
                     p)
  names(info)[4:5] <- c("from", "to")
  
  
  x <- do.call(cbind.DNAbin, c(x, fill.with.gaps = TRUE))
  ntax <- (nrow(x))
  nchar <- sum(info$length)
  #   if ( datatype == "distances" ){
  #     ntax <- attr(x, "Size")
  #     taxnames <- labels(x)
  #   } else {
  #     ntax <- nrow(x)
  
  #   }
  
  ## assemble NEXUS file
  ## -------------------
  nex <- c(
    "#NEXUS", 
    paste("\n[created by ips on ", 
          date(), "]\n", sep = ""))
  
  # TAXA BLOCK (optional)
  # ---------------------
  if ( taxblock ){
    nex <- c(
      nex,
      "begin taxa;", 
      paste("\tdimensions ntax=", ntax,";", sep = ""),
      "\ttaxlabels",
      paste("\t", rownames(x), sep = ""), ";\n"
    )
  }
  
  ## adding whitespace to taxonnames to 
  ## get equal string lengths
  ## ------------------------
  ws <- nchar(rownames(x))
  ws <- max(ws) - ws
  ws <- lapply(ws, function(x) paste(rep(" ", x), 
                                     collapse = ""))
  rownames(x) <- paste(rownames(x), ws, sep = "")
  
  if ( is.numeric(block.width ) ){
    interleave <- " interleave"
  } else {
    interleave <- ifelse(nrow(info) > 1, " interleave", "")
    block.width <- NULL
  }
  
  m <- vector("list", nrow(info))
  for ( i in seq_along(m) ){
    mm <- x[, info$from[i]:info$to[i]]
    if ( is.null(block.width) ){
      bw <- ncol(mm)
    } else {
      bw <- block.width
    }
    m[[i]] <- matrixBlock(mm, bw)
    if ( nrow(info) > 1 ){
      cmt <- paste("[Position ", info$from[i], "-", 
                   info$to[i], ": ", rownames(info)[i], 
                   " (", info$length[i], "bp)]", 
                   sep = "")
      m[[i]] <- c(cmt, m[[i]])
    }
    
    
  }
  m <- c("matrix", unlist(m), ";", "end;")
  
  # write DATA BLOCK
  # -----------------
  #   if ( datatype == "dna" || datatype == "standard" ){
  #     if ( datatype == "standard" && !identical(x[1, 1], round(x[1, 1])) ) 
  #       datatype <- "continuous"
  #     if ( !datatype == "standard" ){
  #       dt <- datatype
  #     } else {
  #       dt <- paste(datatype, " symbols=\"", 
  #                   paste(unique(unlist(x)), collapse = " "), "\"", 		
  #                   sep = "")
  #     }
  
  # assemble DATA BLOCK
  # -------------------
  nex <- c(
    nex,
    paste("begin", ifelse(taxblock, "characters;", "data;")),
    paste("\tdimensions ntax=", ntax, 
          " nchar=", nchar, ";", sep = ""),
    paste("\tformat datatype=", unique(info$datatype), 
          " missing=", unique(info$missing),
          " gap=-", interleave, ";", sep = ""),
    m
  )
  #   }
  
  # write NEXUS file
  # ----------------
  if ( missing(file) ) {
    ## return character vector:
    return(nex) 
  } else {
    if ( file == "" ) {
      ## print onto screen:
      cat(nex, sep = "\n")
    } else {
      ## write to file:
      write(nex, file = file)
    }
  }
}
