
library(ggplot2)
library(griph)
library(dplyr)
library(DESeq2)
source("utility_functions.R")

load("results/raw_norm_DESeq_sfn.RData")

res_griph <- griph_cluster(counts(dds), ClassAssignment = dds$exp_group, 
                           BatchAssignment = factor(dds$batch),
                           use.par = FALSE, plot = FALSE)

g.true <- plotGraph(res_griph, fill.type = "true", mark.type = "predicted")

emb <- tibble(x = g.true$y[,1], y = g.true$y[,2], cond = dds$sample_type, 
              cluster = factor(res_griph$MEMB), btch = dds$batch)

col_pal <- c("#7fbc41", "#7fbc41", "#7fbc41", #Epi
             "#6a51a3", "#6a51a3", #PE
             "#C19A6B", #TE
             "#fe9929", "#fe9929", "#fe9929", "#d73027", "#d73027", #AI
             "#de75ae", #mTeSR1 matrigel
             "#c994c7", "#c994c7", #mTeSR1 laminin
             "#1f78b4", "#1f78b4", "#1f78b4", "#1f78b4", "#1f78b4", #KSR
             "#737373", "#737373", "#737373", "#737373", 
             "#737373", "#737373", "#737373", "#737373", "#737373" #Naive
)

shape_pal <- c(16, 17, 15, #Epi
               16, 17, #PE
               16, #TE
               16, 17, 8, 15, 18, #AI
               16, #mTeSR1 matrigel
               15, 18, #mTeSR1 laminin
               16, 16, 17, 15, 18, #KSR
               16, 17, 17, 15, 
               18, 18, 17, 15, 18 #Naive
)

ggplot(emb, aes(x, y, colour = cond, shape = cond)) + 
  geom_point(size = 2.5) +
  scale_colour_manual(values = col_pal) + 
  scale_shape_manual(values = shape_pal) +
  labs(x = "Dim 1", y = "Dim 2", colour = "", shape = "") + 
  theme_bw(base_size = 15)

save(res_griph, emb, col_pal, shape_pal, file = "results/griph_clustering.RData")