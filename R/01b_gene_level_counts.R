
library(tximport)
library(readr)

# Read in the reference sample data
sams <- read_tsv("all_samples.tsv")

# Path to the quantification files
files <- paste("../hESC_project/", sams$batch, "tx_quantification", 
               sams$sample_name, "quant.sf", sep = "/")
names(files) <- sams$sample_name

# Load the transcript to gene ID map
load("../hESC_project/ref_transcriptome/t2g_light.RData")

# Perform the gene-level quantification
hesc <- tximport(files, type = "salmon", tx2gene = t2g)

save(hesc, file = "results/hesc.RData")
