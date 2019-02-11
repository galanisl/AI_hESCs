
library(dplyr)

# Permutation 1
load("../hESC_project/ref_transcriptome/t2g_complete.RData")
load("../AI_hESCs_p1/results/raw_norm_DESeq_sfn.RData")

cData <- tibble(sample_name = dds$sample_name, cell_line = dds$cell_line, 
                medium = dds$medium, rna_seq = dds$rna_seq, 
                technology = dds$technology, batch = dds$batch, 
                ENA_study_ID = dds$ENA_study_ID, exp_group = dds$exp_group)

rData <- tibble(symbol = rownames(dds)) %>% 
  left_join(t2g, by = "symbol") %>% 
  select(symbol, gene_id)
rData <- rData %>% 
  filter(!duplicated(symbol))

expr_data_p1 <- SummarizedExperiment(
  assays = list(counts = counts(dds), avgTxLength = assays(dds)$avgTxLength,
                normalizedImputedCounts = sf),
  rowData = rData, 
  colData = cData
  )

write.table(counts(dds), file = "data/expr_matrix_raw_p1.tsv", quote = FALSE, sep = "\t")
write.table(sf, file = "data/expr_matrix__norm_p1.tsv", quote = FALSE, sep = "\t")
saveRDS(expr_data_p1, file = "data/expr_data_all_p1.RDS")

rm(list = ls())

# Permutation 2
load("../hESC_project/ref_transcriptome/t2g_complete.RData")
load("../AI_hESCs_p2/results/raw_norm_DESeq_sfn.RData")

cData <- tibble(sample_name = dds$sample_name, cell_line = dds$cell_line, 
                medium = dds$medium, rna_seq = dds$rna_seq, 
                technology = dds$technology, batch = dds$batch, 
                ENA_study_ID = dds$ENA_study_ID, exp_group = dds$exp_group)

rData <- tibble(symbol = rownames(dds)) %>% 
  left_join(t2g, by = "symbol") %>% 
  select(symbol, gene_id)
rData <- rData %>% 
  filter(!duplicated(symbol))

expr_data_p2 <- SummarizedExperiment(
  assays = list(counts = counts(dds), avgTxLength = assays(dds)$avgTxLength,
                normalizedImputedCounts = sf),
  rowData = rData, 
  colData = cData
)

write.table(counts(dds), file = "data/expr_matrix_raw_p2.tsv", quote = FALSE, sep = "\t")
write.table(sf, file = "data/expr_matrix__norm_p2.tsv", quote = FALSE, sep = "\t")
saveRDS(expr_data_p2, file = "data/expr_data_all_p2.RDS")
