
library(DESeq2)
library(dplyr)
library(readr)
source("utility_functions.R")

load("results/hesc.RData")

stype <- read_tsv("all_samples.tsv")
stype <- stype %>% 
  mutate(general_group = parse_factor(general_group, NULL))

# Create a DESeq2 object with a more general design
# Putting all the data together
dds <- DESeqDataSetFromTximport(txi = hesc, colData = stype, 
                                design = ~ general_group)

mito <- grep("^MT-", rownames(dds))
dds <- dds[-mito, ]
ribo <- read_tsv("ref_transcriptome/ribosomal.txt")
dds <- dds[!(rownames(dds) %in% ribo$`Approved Symbol`), ]
load("ref_transcriptome/t2g_complete.RData")
pseudo <- t2g[grep("pseudogene", t2g$gene_biotype), ]
dds <- dds[!(rownames(dds) %in% pseudo$symbol), ]
dds <- dds[!near(rowSums(counts(dds)), 0), ]
dds <- estimateSizeFactors(dds)

# DESeq2 differential expression analysis
dds <- DESeq(dds)
save(dds, file = "results/DESeq_with_contrasts.RData")

# Obtain contrasts of interest and coerce to tibble. The p-values correspond to
# the null hypothesis that the expected log2FCs are >= 0
ai_v_epi <- results(dds, contrast = c("general_group", "AI", "Epi"), lfcThreshold = 1)
ai_v_epi <- DESeqRes2tibble(ai_v_epi)

mtesr_v_epi <- results(dds, contrast = c("general_group", "mTeSR1", "Epi"), lfcThreshold = 1)
mtesr_v_epi <- DESeqRes2tibble(mtesr_v_epi)

ksr_v_epi <- results(dds, contrast = c("general_group", "KSR", "Epi"), lfcThreshold = 1)
ksr_v_epi <- DESeqRes2tibble(ksr_v_epi)

t2il_v_epi <- results(dds, contrast = c("general_group", "t2iL_Go", "Epi"), lfcThreshold = 1)
t2il_v_epi <- DESeqRes2tibble(t2il_v_epi)

# Write results to TSV files
write_tsv(ai_v_epi, path = "results/AI_v_Epi.tsv")
write_tsv(mtesr_v_epi, path = "results/mTeSR1_vs_Epi.tsv")
write_tsv(ksr_v_epi, path = "results/KSR_v_Epi.tsv")
write_tsv(t2il_v_epi, path = "results/t2iL_v_Epi.tsv")

# Identify up- and down-regulated genes based on stringent thresholds and 
# perform functional enrichment analysis
library(FunEnrich)
library(cowplot)
pval <- 0.001
log2FC <- log2(1)

updw_ai <- funenrich_analysis(ai_v_epi, log2FC, pval, benj = FALSE)
updw_ksr <- funenrich_analysis(ksr_v_epi, log2FC, pval, benj = FALSE)
updw_mtesr <- funenrich_analysis(mtesr_v_epi, log2FC, pval, benj = FALSE)
updw_t2il <- funenrich_analysis(t2il_v_epi, log2FC, pval, benj = FALSE)

save(updw_ai, updw_ksr, updw_mtesr, updw_t2il, file = "results/perm2_DGE.RData")

# AI
readr::write_tsv(updw_ai$enr_up$bp, path = "results/up_ai_bp.tsv")
readr::write_tsv(updw_ai$enr_up$cc, path = "results/up_ai_cc.tsv")
readr::write_tsv(updw_ai$enr_up$mf, path = "results/up_ai_mf.tsv")
readr::write_tsv(updw_ai$enr_up$reactome, path = "results/up_ai_reac.tsv")

readr::write_tsv(updw_ai$enr_dw$bp, path = "results/dw_ai_bp.tsv")
readr::write_tsv(updw_ai$enr_dw$cc, path = "results/dw_ai_cc.tsv")
readr::write_tsv(updw_ai$enr_dw$mf, path = "results/dw_ai_mf.tsv")
readr::write_tsv(updw_ai$enr_dw$reactome, path = "results/dw_ai_reac.tsv")

# KSR and 3iL
readr::write_tsv(updw_ksr$enr_up$bp, path = "results/up_ksr_bp.tsv")
readr::write_tsv(updw_ksr$enr_up$cc, path = "results/up_ksr_cc.tsv")
readr::write_tsv(updw_ksr$enr_up$mf, path = "results/up_ksr_mf.tsv")
readr::write_tsv(updw_ksr$enr_up$reactome, path = "results/up_ksr_reac.tsv")

readr::write_tsv(updw_ksr$enr_dw$bp, path = "results/dw_ksr_bp.tsv")
readr::write_tsv(updw_ksr$enr_dw$cc, path = "results/dw_ksr_cc.tsv")
readr::write_tsv(updw_ksr$enr_dw$mf, path = "results/dw_ksr_mf.tsv")
readr::write_tsv(updw_ksr$enr_dw$reactome, path = "results/dw_ksr_reac.tsv")

# mTeSR1
readr::write_tsv(updw_mtesr$enr_up$bp, path = "results/up_mtesr_bp.tsv")
readr::write_tsv(updw_mtesr$enr_up$cc, path = "results/up_mtesr_cc.tsv")
readr::write_tsv(updw_mtesr$enr_up$mf, path = "results/up_mtesr_mf.tsv")
readr::write_tsv(updw_mtesr$enr_up$reactome, path = "results/up_mtesr_reac.tsv")

readr::write_tsv(updw_mtesr$enr_dw$bp, path = "results/dw_mtesr_bp.tsv")
readr::write_tsv(updw_mtesr$enr_dw$cc, path = "results/dw_mtesr_cc.tsv")
readr::write_tsv(updw_mtesr$enr_dw$mf, path = "results/dw_mtesr_mf.tsv")
readr::write_tsv(updw_mtesr$enr_dw$reactome, path = "results/dw_mtesr_reac.tsv")

# t2iL+GO and 5iLA
readr::write_tsv(updw_t2il$enr_up$bp, path = "results/up_t2il_bp.tsv")
readr::write_tsv(updw_t2il$enr_up$cc, path = "results/up_t2il_cc.tsv")
readr::write_tsv(updw_t2il$enr_up$mf, path = "results/up_t2il_mf.tsv")
readr::write_tsv(updw_t2il$enr_up$reactome, path = "results/up_t2il_reac.tsv")

readr::write_tsv(updw_t2il$enr_dw$bp, path = "results/dw_t2il_bp.tsv")
readr::write_tsv(updw_t2il$enr_dw$cc, path = "results/dw_t2il_cc.tsv")
readr::write_tsv(updw_t2il$enr_dw$mf, path = "results/dw_t2il_mf.tsv")
readr::write_tsv(updw_t2il$enr_dw$reactome, path = "results/dw_t2il_reac.tsv")




