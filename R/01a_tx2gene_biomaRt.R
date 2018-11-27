
# This scripts converts transcript IDs to Ensembl gene IDs and official gene
# symbol with the biomaRt tool

library(biomaRt)
library(dplyr)

# Read one of Salmon's quantification files to obtain the list of transcript IDs
tx <- readr::read_tsv("wamaitha/tx_quantification/WAM464A1_S1/quant.sf") %>% 
  select(Name) %>% 
  rename(tx_idv = Name)

# Remove the version number from the transcript IDs
tx <- tx %>% 
  mutate(tx_id = sapply(strsplit(tx_idv, "[.]"), `[`, 1))

# Map the transcript IDs to gene IDs
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")

g_ids <- getBM(attributes = c("ensembl_transcript_id", "hgnc_symbol", 
                              "ensembl_gene_id", "gene_biotype"), 
               filters = "ensembl_transcript_id", 
               values = tx$tx_id, 
               mart = ensembl)

g_ids <- as_tibble(g_ids) %>% 
  rename(tx_id = ensembl_transcript_id, symbol = hgnc_symbol, 
         gene_id = ensembl_gene_id)

# Join the transcript and gene tables to form the final map
t2g <- left_join(tx, g_ids, by = "tx_id") %>% 
  filter(!duplicated(tx_id, fromLast = TRUE))

# Substitute unmapped (NA) cases with their transcript ID
t2g$symbol[is.na(t2g$symbol)] <- t2g$tx_id[is.na(t2g$symbol)]
t2g$gene_id[is.na(t2g$gene_id)] <- t2g$tx_id[is.na(t2g$gene_id)]
t2g$gene_biotype[is.na(t2g$gene_biotype)] <- "unknown"

# Substitute unmapped ('') cases with their gene ID
t2g$symbol[t2g$symbol == ""] <- t2g$gene_id[t2g$symbol == ""]

# Save the complete map and a light version with the mapping of interest
save(t2g, file = "ref_transcriptome/t2g_complete.RData")

t2g <- t2g %>% 
  select(tx_idv, symbol)

save(t2g, file = "ref_transcriptome/t2g_light.RData")

