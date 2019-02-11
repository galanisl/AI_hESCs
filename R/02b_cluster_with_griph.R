
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
              cluster = factor(res_griph$MEMB))

<<<<<<< HEAD
# col_pal <- c("#4d9221", "#7fbc41", "#a6d96a", #Epi
#              "#9e9ac8", "#6a51a3", #PE
#              "#C19A6B", #TE
#              "#a50026", "#d73027", "#f46d43", "#fdae61", #AI
#              "#f1b6da", "#de77ae", #mTeSR1
#              "#1f78b4", "#6baed6", "#41b6c4", "#1d91c0", "#225ea8", "#7fcdbb", #KSR
#              "#d9d9d9", "#bdbdbd", "#969696", "#737373", 
#              "#525252", "#464646", "#252525", "#101010", "#000000" #Naive
# )

col_pal <- c("#7fbc41", "#7fbc41", "#7fbc41", #Epi
             "#6a51a3", "#6a51a3", #PE
             "#C19A6B", #TE
             "#d73027", "#d73027", "#d73027", "#d73027", #AI
             "#de75ae", "#de75ae", #mTeSR1
             "#1f78b4", "#1f78b4", "#1f78b4", "#1f78b4", "#1f78b4", "#9ecae1", #KSR
             "#737373", "#737373", "#737373", "#737373", 
             "#737373", "#737373", "#737373", "#737373", "#737373" #Naive
)

shape_pal <- c(16, 17, 15, #Epi
               16, 17, #PE
               16, #TE
               16, 17, 15, 18, #AI
               16, 17, #mTeSR1
               16, 16, 17, 15, 18, 16, #KSR
               16, 17, 17, 15, 
               18, 18, 17, 15, 18 #Naive
)

ggplot(emb, aes(x, y, colour = cond, shape = cond)) + 
  geom_point(size = 2.5) +
  scale_colour_manual(values = col_pal) + 
  scale_shape_manual(values = shape_pal) +
=======
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
>>>>>>> e1a244ea72e6fc0478d258ff4d5ff9198c8dc9b2
  labs(x = "Dim 1", y = "Dim 2", colour = "") + 
  theme_bw(base_size = 15)

save(res_griph, emb, col_pal, file = "results/griph_clustering.RData")
