#' @title PATHd8
#' @description This function is a wrapper for PATHd8 and can be used for
#'   phylogenetic dating, especially of large trees
#' @param phy An object of class \code{\link[ape:read.tree]{phylo}}.
#' @param exec A character string giving the path to the PATHd8 program.
#' @param seql sequence length of alignment
#' @param calibration A data frame with 4 columns and as many rows as 
#'   calibration points. Columns are: taxon 1; taxon 2; one of c("minage", 
#'   "maxage", "fixage"); age.
#' @return tree list of ultrametric trees returned from PATHd8 of class 
#'   \code{\link[ape:read.tree]{phylo}}. First tree is PATHd8 chronogram, which is a calibrated
#'   ultrametric tree. Second is a PATH tree, which is a ultrametric tree 
#'   without calibration.
#' @references Britton et al (2006). PATHd8---a new method for estimating 
#'   divergence times in large phylogenetic trees without a molecular clock. 
#'   Available from the authors (www.math.su.se/PATHd8)
#' @references Britton et al. (2007). Estimating divergence times in large 
#'   phylogenetic trees. Systematic biology. 56:741--752
#' @author Franz-Sebastian Krah
#' @examples \dontrun{
#' ## This example is taken from the PATHD8 manual
#' cal <- rbind(cal1 = c("Rat", "Ostrich", "minage", 260), 
#' cal2 = c("Human", "Platypus", "fixage", 125),
#' cal3 = c("Alligator", "Ostrich", "minage", 150))
#' colnames(cal) = c("tax1", "tax2", "age_type", "age")
#' phy <- read.tree(text = paste0("((((Rat:0.007148,Human:0.001808):0.024345,",
#'                                "Platypus:0.016588):0.012920,(Ostrich:0.018119,",
#'                                "Alligator:0. 006232):0.004708):0.028037,Frog:0);")
#' seql <- 1823
#' pathd8(phy, exec = "/Applications/PATHd8/PATHd8", seql = seql, calibration = cal)
#' }
#' @importFrom plyr rbind.fill
#' @export


pathd8 <- function(phy, exec = "/Applications/PATHd8/PATHd8", seql, calibration){
  
  if (!inherits(phy, "phylo"))
    stop("'phy' is not of class phylo")
  
  fns <- vector(length = 6)
  for (i in seq_along(fns))
    fns[i] <- tempfile(pattern = "pathd8", tmpdir = tempdir(), fileext = ".nwk")
  unlink(fns[file.exists(fns)])

  tree <- write.tree(phy)
  cal <- calibration 
  cal <- paste0("mrca: ", cal[,1], ", ", cal[,2], ", ", cal[,3], "=", cal[,4], ";")
  
  if (!missing(calibration)){
    write(paste0("Sequence length= ", 
                 seql, ";", 
                 " \n \n", 
                 tree, 
                 "\n \n",
                 paste(cal, collapse = "\n")),
          file = fns[1])
  }
  # write.tree(phy, file = fns[1])
  readLines(fns[1])
  call <- paste(exec, fns[1], fns[2])
  system(call)
  
  # pathd8 tree
  out <- readLines(fns[2])
  trst1 <- grep("d8 tree", out)
  tr1 <- out[trst1]
  write(tr1, file = fns[3])
  tr <- read.tree(fns[3])
  
  ## MPL tree
  trst2 <- grep("MPL tree", out)
  tr2 <- out[trst2]
  write(tr2, file = fns[4])
  tr2 <- read.tree(fns[4])
  
  ## Tables
  st_tab1 <- trst1 + 2
  end_tab1 <- st_tab1 + tr$Nnode
  tab1 <- gsub("[[:space:]]+", " ", out[(st_tab1 + 1):end_tab1])
  tab1 <- do.call(rbind, strsplit(tab1, " "))
  tab1 <- data.frame(tab1)
  header <- out[st_tab1]
  header <- gsub("\\b\\s\\b", "_", header)
  header <- gsub("[[:space:]]+", " ", header)
  header <- gsub("\\*\\s", "", header)
  colnames(tab1) <- strsplit(header, " ")[[1]]
  
  ## Table MPL
  st_tab1 <- trst2 + 2
  end_tab1 <- st_tab1 + tr$Nnode
  tab2 <- gsub("[[:space:]]+", " ", out[(st_tab1 + 1):end_tab1])
  tab2 <- rbind.fill(lapply(strsplit(tab2, " "), function(y) { as.data.frame(t(y), stringsAsFactors=FALSE) }))
  tab2 <- data.frame(tab2)
  header <- out[st_tab1]
  header <- gsub("\\b\\s\\b", "_", header)
  header <- gsub("[[:space:]]+", " ", header)
  header <- gsub("\\*\\s", "", header)
  header <- strsplit(header, " ")[[1]]
  # header <- header[-(6:7)]
  colnames(tab2) <- c(header[1:4], c("plumin", "MPLvar"), header[5], c("test", "sig"))
  # tab2 <- tab2[, -8]
  
  ## Delete temporary files
  files <- list.files(tempdir(), full.names = TRUE)
  files <- files[-grep("rs-graphics", files)]
  unlink(files, force = TRUE, recursive = TRUE)
  
  return(list(path_tree = tr, mpl_tree = tr2, path_tab = tab1, mpl_tab = tab2))
}
