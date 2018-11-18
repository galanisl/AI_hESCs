
library(ggplot2)
library(griph)
library(dplyr)
library(DESeq2)
source("utility_functions.R")

load("results/raw_norm_DESeq_sfn.RData")

# dds <- dds[, colData(dds)$condition %in% c("Epi", "TE", "PE")]
# dds <- dds[, colData(dds)$rna_seq == "single_cell"]

# Determine the fraction of variable genes to retain based on robust z-scores
# z-score = 1.65 -> p-val = 0.05
# meanDM <- mean(counts(dds))
# nSD <- apply(counts(dds), 1, function(x) sd(x)/meanDM)
# ConstRows <- which(nSD < 1e-3)
# 
# zscores <- select_variable_genes(counts(dds)[-ConstRows, ])
# frac <- sum(zscores >= 1.65)/nrow(counts(dds)[-ConstRows, ])

res_griph <- griph_cluster(counts(dds), ClassAssignment = dds$exp_group, 
                     BatchAssignment = factor(dds$batch),
                     use.par = FALSE, plot = FALSE)

g.true <- plotGraph(res_griph, fill.type = "true", mark.type = "predicted")

emb <- tibble(x = g.true$y[,1], y = g.true$y[,2], cond = dds$sample_type, 
              cluster = factor(res_griph$MEMB))

col_pal <- c("#4d9221", "#7fbc41", "#a6d96a", #Epi
             "#9e9ac8", "#6a51a3", #PE
             "#C19A6B", #TE
             "#a50026", "#d73027", "#f46d43", "#fdae61", #AI
             "#f1b6da", "#de77ae", #mTeSR1
             "#1f78b4", "#6baed6", "#41b6c4", "#1d91c0", "#225ea8", "#7fcdbb", #KSR
             "#d9d9d9", "#bdbdbd", "#969696", "#737373", 
             "#525252", "#464646", "#252525", "#101010", "#000000" #Naive
)

ggplot(emb, aes(x, y, colour = cond)) + 
  geom_point(size = 2.5) +
  scale_colour_manual(values = col_pal) + 
  labs(x = "Dim 1", y = "Dim 2", colour = "") + 
  theme_bw(base_size = 15)

save(res_griph, emb, col_pal, file = "results/griph_clustering.RData")
