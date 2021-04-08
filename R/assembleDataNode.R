## This code is part of the ips package
## © C. Heibl 2014 (last update 2020-03-13)

#' @export

assembleDataNode <- function(DNAbin){
  
  data <- vector(mode = "list")
  for (i in 1:length(DNAbin)){
    ss  <-  as.list(DNAbin[[i]])
    ss <- lapply(as.character(ss), paste, collapse = "")
    ss <- lapply(ss, gsub, pattern = "n", replacement = "?")
    id <- paste("seq", names(ss), sep = "_")
    if (i > 1) {
      id <- paste(id, i, sep = "_")
    }
    ss <- data.frame(id = id,
                     spec = "Sequence",
                     taxon = names(ss),
                     totalcount = "4", ## DNA
                     value = toupper(unlist(ss)),
                     stringsAsFactors = FALSE)
    id <- names(DNAbin)[i]
    # id <- paste("BEAST", id, sep = "_") has to be adapted throughout XML
    sequence <- apply(ss, 1, function(x) xmlNode("sequence", 
                                                 attrs = c(x[1], x[2], x[3], x[4], x[5])))
    attr <- c(id = id, spec = "Alignment")
    if (i == 1) attr <- c(attr, name = "alignment")
    data <- c(data, list(xmlNode("data", attrs = attr,
                                 .children = sequence)))
  }
  data
}
