
library(tximport)
library(readr)

# Read in the reference sample data (this file is in folder 'report' on GitHub)
sams <- read_tsv("samples_features.tsv")

# Path to the quantification files
files <- paste(sams$batch, "tx_quantification", sams$sample_name, "quant.sf", 
               sep = "/")
names(files) <- sams$sample_name

# Load the transcript to gene ID map
load("ref_transcriptome/t2g_light.RData")

# Perform the gene-level quantification
hesc <- tximport(files, type = "salmon", tx2gene = t2g)

save(hesc, file = "results/hesc.RData")
