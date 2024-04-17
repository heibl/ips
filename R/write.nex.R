## This code is part of the ips package
## © C. Heibl 2014 (last update 2019-06-21)

#' @export

write.nex <- function(x, file, block.width = 60, 
                      taxblock = FALSE){
  
  ## A data frame is a list: is.list(data.frame) == TRUE !!!
  if (!is.list(x) | is.data.frame(x)){
    x <- list(x)
  }
    
  ## Auxiliary function 1: Asses datatype
  ## ------------------------------------
  getDataType <- function(x){
    datatype <- class(x)
    datatype[datatype == "DNAbin"] <- "dna"
    datatype[datatype == "dist"] <- "distances"
    datatype[datatype == "data.frame"] <- "standard"
    datatype
  }
  

  ## Auxiliary function 2: Asses token used for missing data 
  ## (function adapted for data frames 2016-01-26)
  ## ---------------------------------------------
  m <- function(x, datatype) {
    if (getDataType(x) == "standard") {
      n <- ifelse(any(x == "?"), "?", "N")
    } else {
      n <- ifelse("as.raw(2)" %in% x, "?", "N")
    }
    n
  }
  ## Nucleotide positions in partitions
  ## ----------------------------------
  p <- cbind(rep(1, length(x)), sapply(x, ncol))
  if (nrow(p) > 1) {
    for ( i in 2:nrow(p) ){
      p[i, 1] <- p[i - 1, 2] + 1
      p[i, 2] <- p[i, 1] + p[i, 2] -1
    }
  }
  
  ## Assemble info data
  ## ------------------
  info <- data.frame(datatype = sapply(x, getDataType),
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
  
  if ( is.numeric(block.width) ){
    interleave <- " interleave"
  } else {
    interleave <- ifelse(nrow(info) > 1, " interleave", "")
    block.width <- NULL
  }
  
  m <- vector("list", nrow(info))
  for (i in seq_along(m)){
    mm <- x[, info$from[i]:info$to[i]]
    bw <- ifelse(is.null(block.width), ncol(mm), block.width)
    m[[i]] <- matrixBlock(mm, bw)
    if (nrow(info) > 1){
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
