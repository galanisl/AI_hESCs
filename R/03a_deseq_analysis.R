
library(DESeq2)
library(dplyr)
library(readr)
source("utility_functions.R")

<<<<<<< HEAD
load("results/hesc_all.RData")
=======
load("results/hesc.RData")
>>>>>>> 4cab0b77e62cced71c7460c1a0d9f17c9751de04

stype <- read_tsv("all_samples.tsv")
stype <- stype %>% 
  mutate(general_group = parse_factor(general_group, NULL))

# Create a DESeq2 object with a more general design
# Putting all the data together
<<<<<<< HEAD
hesc$abundance <- hesc$abundance[, stype$sample_name]
hesc$counts <- hesc$counts[, stype$sample_name]
hesc$length <- hesc$length[, stype$sample_name]
=======
>>>>>>> 4cab0b77e62cced71c7460c1a0d9f17c9751de04
dds <- DESeqDataSetFromTximport(txi = hesc, colData = stype, 
                                design = ~ general_group)

mito <- grep("^MT-", rownames(dds))
dds <- dds[-mito, ]
<<<<<<< HEAD
ribo <- read_tsv("data/ribosomal.txt")
dds <- dds[!(rownames(dds) %in% ribo$`Approved Symbol`), ]
load("data/t2g_complete.RData")
=======
ribo <- read_tsv("ref_transcriptome/ribosomal.txt")
dds <- dds[!(rownames(dds) %in% ribo$`Approved Symbol`), ]
load("ref_transcriptome/t2g_complete.RData")
>>>>>>> 4cab0b77e62cced71c7460c1a0d9f17c9751de04
pseudo <- t2g[grep("pseudogene", t2g$gene_biotype), ]
dds <- dds[!(rownames(dds) %in% pseudo$symbol), ]
dds <- dds[!near(rowSums(counts(dds)), 0), ]
dds <- estimateSizeFactors(dds)

# DESeq2 differential expression analysis
dds <- DESeq(dds)
save(dds, file = "results/DESeq_with_contrasts.RData")

# Obtain contrasts of interest and coerce to tibble. The p-values correspond to
<<<<<<< HEAD
# the null hypothesis that the expected log2FCs are = 0
# Epi
ai_v_epi <- results(dds, contrast = c("general_group", "AI", "Epi"), lfcThreshold = 0)
ai_v_epi <- DESeqRes2tibble(ai_v_epi)

mtesr_v_epi <- results(dds, contrast = c("general_group", "mTeSR1", "Epi"), lfcThreshold = 0)
mtesr_v_epi <- DESeqRes2tibble(mtesr_v_epi)

mtesrL_v_epi <- results(dds, contrast = c("general_group", "mTeSR1_laminin", "Epi"), lfcThreshold = 0)
mtesrL_v_epi <- DESeqRes2tibble(mtesrL_v_epi)

ksr_v_epi <- results(dds, contrast = c("general_group", "KSR", "Epi"), lfcThreshold = 0)
ksr_v_epi <- DESeqRes2tibble(ksr_v_epi)

t2il_v_epi <- results(dds, contrast = c("general_group", "t2iL_Go", "Epi"), lfcThreshold = 0)
t2il_v_epi <- DESeqRes2tibble(t2il_v_epi)

# TE
ai_v_te <- results(dds, contrast = c("general_group", "AI", "TE"), lfcThreshold = 0)
ai_v_te <- DESeqRes2tibble(ai_v_te)

mtesr_v_te <- results(dds, contrast = c("general_group", "mTeSR1", "TE"), lfcThreshold = 0)
mtesr_v_te <- DESeqRes2tibble(mtesr_v_te)

mtesrL_v_te <- results(dds, contrast = c("general_group", "mTeSR1_laminin", "TE"), lfcThreshold = 0)
mtesrL_v_te <- DESeqRes2tibble(mtesrL_v_te)

ksr_v_te <- results(dds, contrast = c("general_group", "KSR", "TE"), lfcThreshold = 0)
ksr_v_te <- DESeqRes2tibble(ksr_v_te)

t2il_v_te <- results(dds, contrast = c("general_group", "t2iL_Go", "TE"), lfcThreshold = 0)
t2il_v_te <- DESeqRes2tibble(t2il_v_te)

#PE
ai_v_pe <- results(dds, contrast = c("general_group", "AI", "PE"), lfcThreshold = 0)
ai_v_pe <- DESeqRes2tibble(ai_v_pe)

mtesr_v_pe <- results(dds, contrast = c("general_group", "mTeSR1", "PE"), lfcThreshold = 0)
mtesr_v_pe <- DESeqRes2tibble(mtesr_v_pe)

mtesrL_v_pe <- results(dds, contrast = c("general_group", "mTeSR1_laminin", "PE"), lfcThreshold = 0)
mtesrL_v_pe <- DESeqRes2tibble(mtesrL_v_pe)

ksr_v_pe <- results(dds, contrast = c("general_group", "KSR", "PE"), lfcThreshold = 0)
ksr_v_pe <- DESeqRes2tibble(ksr_v_pe)

t2il_v_pe <- results(dds, contrast = c("general_group", "t2iL_Go", "PE"), lfcThreshold = 0)
t2il_v_pe <- DESeqRes2tibble(t2il_v_pe)

# Blastocyst
epi_v_pe <- results(dds, contrast = c("general_group", "Epi", "PE"), lfcThreshold = 0)
epi_v_pe <- DESeqRes2tibble(epi_v_pe)

epi_v_te <- results(dds, contrast = c("general_group", "Epi", "TE"), lfcThreshold = 0)
epi_v_te <- DESeqRes2tibble(epi_v_te)

te_v_pe <- results(dds, contrast = c("general_group", "TE", "PE"), lfcThreshold = 0)
te_v_pe <- DESeqRes2tibble(te_v_pe)

# Write results to TSV files
write_tsv(ai_v_epi, path = "dgea/AI_v_Epi.tsv")
write_tsv(mtesr_v_epi, path = "dgea/mTeSR1_vs_Epi.tsv")
write_tsv(mtesrL_v_epi, path = "dgea/mTeSR1L_vs_Epi.tsv")
write_tsv(ksr_v_epi, path = "dgea/KSR_v_Epi.tsv")
write_tsv(t2il_v_epi, path = "dgea/t2iL_v_Epi.tsv")

write_tsv(ai_v_te, path = "dgea/AI_v_TE.tsv")
write_tsv(mtesr_v_te, path = "dgea/mTeSR1_vs_TE.tsv")
write_tsv(mtesrL_v_te, path = "dgea/mTeSR1L_vs_TE.tsv")
write_tsv(ksr_v_te, path = "dgea/KSR_v_TE.tsv")
write_tsv(t2il_v_te, path = "dgea/t2iL_v_TE.tsv")

write_tsv(ai_v_pe, path = "dgea/AI_v_PE.tsv")
write_tsv(mtesr_v_pe, path = "dgea/mTeSR1_vs_PE.tsv")
write_tsv(mtesrL_v_pe, path = "dgea/mTeSR1L_vs_PE.tsv")
write_tsv(ksr_v_pe, path = "dgea/KSR_v_PE.tsv")
write_tsv(t2il_v_pe, path = "dgea/t2iL_v_PE.tsv")

write_tsv(epi_v_pe, path = "dgea/Epi_v_PE.tsv")
write_tsv(epi_v_te, path = "dgea/Epi_v_TE.tsv")
write_tsv(te_v_pe, path = "dgea/TE_v_PE.tsv")
=======
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
>>>>>>> 4cab0b77e62cced71c7460c1a0d9f17c9751de04

# Identify up- and down-regulated genes based on stringent thresholds and 
# perform functional enrichment analysis
library(FunEnrich)
library(cowplot)
<<<<<<< HEAD
library(ggplot2)
pval <- 0.001
log2FC <- log2(2)
=======
pval <- 0.001
log2FC <- log2(1)
>>>>>>> 4cab0b77e62cced71c7460c1a0d9f17c9751de04

updw_ai <- funenrich_analysis(ai_v_epi, log2FC, pval, benj = FALSE)
updw_ksr <- funenrich_analysis(ksr_v_epi, log2FC, pval, benj = FALSE)
updw_mtesr <- funenrich_analysis(mtesr_v_epi, log2FC, pval, benj = FALSE)
<<<<<<< HEAD
updw_mtesrL <- funenrich_analysis(mtesrL_v_epi, log2FC, pval, benj = FALSE)
updw_t2il <- funenrich_analysis(t2il_v_epi, log2FC, pval, benj = FALSE)

save(updw_ai, updw_ksr, updw_mtesr, updw_mtesrL, updw_t2il, file = "dgea/dgea.RData")

# AI
readr::write_tsv(updw_ai$enr_up$bp, path = "dgea/up_ai_bp.tsv")
readr::write_tsv(updw_ai$enr_up$cc, path = "dgea/up_ai_cc.tsv")
readr::write_tsv(updw_ai$enr_up$mf, path = "dgea/up_ai_mf.tsv")
readr::write_tsv(updw_ai$enr_up$reactome, path = "dgea/up_ai_reac.tsv")

readr::write_tsv(updw_ai$enr_dw$bp, path = "dgea/dw_ai_bp.tsv")
readr::write_tsv(updw_ai$enr_dw$cc, path = "dgea/dw_ai_cc.tsv")
readr::write_tsv(updw_ai$enr_dw$mf, path = "dgea/dw_ai_mf.tsv")
readr::write_tsv(updw_ai$enr_dw$reactome, path = "dgea/dw_ai_reac.tsv")

# KSR/FGF
readr::write_tsv(updw_ksr$enr_up$bp, path = "dgea/up_ksr_bp.tsv")
readr::write_tsv(updw_ksr$enr_up$cc, path = "dgea/up_ksr_cc.tsv")
readr::write_tsv(updw_ksr$enr_up$mf, path = "dgea/up_ksr_mf.tsv")
readr::write_tsv(updw_ksr$enr_up$reactome, path = "dgea/up_ksr_reac.tsv")

readr::write_tsv(updw_ksr$enr_dw$bp, path = "dgea/dw_ksr_bp.tsv")
readr::write_tsv(updw_ksr$enr_dw$cc, path = "dgea/dw_ksr_cc.tsv")
readr::write_tsv(updw_ksr$enr_dw$mf, path = "dgea/dw_ksr_mf.tsv")
readr::write_tsv(updw_ksr$enr_dw$reactome, path = "dgea/dw_ksr_reac.tsv")

# mTeSR1 matrigel
readr::write_tsv(updw_mtesr$enr_up$bp, path = "dgea/up_mtesr_bp.tsv")
readr::write_tsv(updw_mtesr$enr_up$cc, path = "dgea/up_mtesr_cc.tsv")
readr::write_tsv(updw_mtesr$enr_up$mf, path = "dgea/up_mtesr_mf.tsv")
readr::write_tsv(updw_mtesr$enr_up$reactome, path = "dgea/up_mtesr_reac.tsv")

readr::write_tsv(updw_mtesr$enr_dw$bp, path = "dgea/dw_mtesr_bp.tsv")
readr::write_tsv(updw_mtesr$enr_dw$cc, path = "dgea/dw_mtesr_cc.tsv")
readr::write_tsv(updw_mtesr$enr_dw$mf, path = "dgea/dw_mtesr_mf.tsv")
readr::write_tsv(updw_mtesr$enr_dw$reactome, path = "dgea/dw_mtesr_reac.tsv")

# mTeSR1 laminin
readr::write_tsv(updw_mtesrL$enr_up$bp, path = "dgea/up_mtesrL_bp.tsv")
readr::write_tsv(updw_mtesrL$enr_up$cc, path = "dgea/up_mtesrL_cc.tsv")
readr::write_tsv(updw_mtesrL$enr_up$mf, path = "dgea/up_mtesrL_mf.tsv")
readr::write_tsv(updw_mtesrL$enr_up$reactome, path = "dgea/up_mtesrL_reac.tsv")

readr::write_tsv(updw_mtesrL$enr_dw$bp, path = "dgea/dw_mtesrL_bp.tsv")
readr::write_tsv(updw_mtesrL$enr_dw$cc, path = "dgea/dw_mtesrL_cc.tsv")
readr::write_tsv(updw_mtesrL$enr_dw$mf, path = "dgea/dw_mtesrL_mf.tsv")
readr::write_tsv(updw_mtesrL$enr_dw$reactome, path = "dgea/dw_mtesrL_reac.tsv")

# t2iL+GO and 5iLA
readr::write_tsv(updw_t2il$enr_up$bp, path = "dgea/up_t2il_bp.tsv")
readr::write_tsv(updw_t2il$enr_up$cc, path = "dgea/up_t2il_cc.tsv")
readr::write_tsv(updw_t2il$enr_up$mf, path = "dgea/up_t2il_mf.tsv")
readr::write_tsv(updw_t2il$enr_up$reactome, path = "dgea/up_t2il_reac.tsv")

readr::write_tsv(updw_t2il$enr_dw$bp, path = "dgea/dw_t2il_bp.tsv")
readr::write_tsv(updw_t2il$enr_dw$cc, path = "dgea/dw_t2il_cc.tsv")
readr::write_tsv(updw_t2il$enr_dw$mf, path = "dgea/dw_t2il_mf.tsv")
readr::write_tsv(updw_t2il$enr_dw$reactome, path = "dgea/dw_t2il_reac.tsv")
=======
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




>>>>>>> 4cab0b77e62cced71c7460c1a0d9f17c9751de04
