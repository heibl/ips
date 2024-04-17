## This code is part of the ips package
## © C. Heibl 2019 (last update 2020-04-30)

#' @title Nucleotide-Nucleotide BLAST
#' @description Provides an interface to BLASTN
#' @param query An object of class \code{DNAbin} containing the sequences which
#'   will be blasted against \code{db}.
#' @param db An object of class \code{DNAbin} containing the reference
#'   sequences, i.e. the sequences against which \code{query} will be blasted.
#' @importFrom data.table fread
#' @export

blastn <- function(query, db){
  
  ## Temporary input/output files
  ## ----------------------------
  fns <- c("query", "db", "out")
  fns <- tempfile(pattern = fns, tmpdir = tempdir(), 
                  fileext = c(".fas", ".fas", ".txt"))
  unlink(fns[file.exists(fns)])

  ## Create BLAST database
  ## ---------------------
  write.FASTA(query, fns[2])
  cmd <- paste("/usr/local/ncbi/blast/bin/makeblastdb", 
               "-in", fns[2], 
               "-dbtype nucl",
               "-parse_seqids")
  system(cmd)
  
  ## Prepare Sequences
  ## ------------------
  ## move upstream
  # seqs <- dbReadDNA(x, acc.tab)
  # names(seqs) <- gsub("-", "__", names(seqs)) # hyphen illegal in BLAST!
  write.fas(query, fns[1])
  
  ## Do the BLAST
  ## ------------
  cls <- c("qseqid", "sseqid", 
           "length", 
           "mismatch", 
           "qstart", "qend", "sstart", "send",
           # "gapopen", # Number of gap openings
           "qcovs", # Query Coverage Per Subject
           "pident", # Percentage of identical matches
           "evalue", 
           "bitscore", 
           "sstrand")
  outfmt <- paste0("-outfmt '", paste(c(6, cls), collapse = " "), "'")
  cmd <- paste("/usr/local/ncbi/blast/bin/blastn",
               "-db", fns[2],
               "-query", fns[1],
               "-task blastn", ## traditional BLASTN requiring an exact match of 11
               "-evalue 1000",
               "-max_hsps 1", ## maximum number of HSPs per subject sequence to save for each query
               "-max_target_seqs 10000",
               outfmt,
               "-out", fns[3])
  system(cmd)
  
  ## Parse output 
  ## ------------
  res <- fread(fns[3])
  names(res) <- cls
  
  # ## For every sequence select hit with lowest E-value
  # res <- by(res, res$qseqid, function(z) z[which.min(z$evalue),])
  # res <- do.call(rbind, res)
  
  # res$qseqid <- gsub("__", "-", res$qseqid) ## bring hyphen back!
  # res <- cbind(splitGiTaxon(res$qseqid), res)
  
  unlink(fns)
  res
}