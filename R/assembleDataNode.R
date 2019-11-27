## This code is part of the ips package
## © C. Heibl 2014 (last update 2015-04-05)

assembleDataNode <- function(DNAbin){
  
  data <- vector(mode = "list")
  for (i in 1:length(DNAbin)){
    ss  <-  as.list(DNAbin[[i]])
    ss <- lapply(as.character(ss), paste, collapse = "")
    ss <- lapply(ss, gsub, pattern = "n", replacement = "?")
    ss <- data.frame(id = paste("seq", names(ss), i, sep = "_"),
                     taxon = names(ss),
                     totalcount = "4", ## DNA
                     value = unlist(ss),
                     stringsAsFactors = FALSE)
    id <- names(DNAbin)[i]
    sequence <- apply(ss, 1, function(x) xmlNode("sequence", 
                                                 attrs = c(x[1], x[2], x[3], x[4])))
    data <- c(data, list(xmlNode("data", 
            attrs = c(id = id, name = "alignment"),
            .children = sequence)))
  }
  data
}