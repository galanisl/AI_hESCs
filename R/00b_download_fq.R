
# This script must be run from the bash terminal as follows:
#   Rscript 00b_download_fq.R sample_data.tsv fastq_raw
# The first argument must be a tab-separated file and contain a column called
# fastq_ftp with the FTP address of the samples' FastQ files. The second 
# argument is the name of the folder where these files will be downloaded.

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