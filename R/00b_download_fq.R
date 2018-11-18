
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2){
  stop(paste0("This script requires two arguments: a TSV file with links to ",
              "sample FASTQ files and the name of an output directory"),
       call. = FALSE)
}

# Read the sample data
sdata <- readr::read_tsv(args[1])

# Download the FASTQ files
sapply(1:nrow(sdata), 
       function(x){
         download.file(sdata$fastq_ftp[x], 
                       paste0(args[2], "/", sdata$run_accession[x], 
                              ".fq.gz"))
         }
       )