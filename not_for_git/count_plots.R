
library(ggplot2)
library(scales)
library(dplyr)

load("../AI_hESCs_p2/results/raw_norm_DESeq_sfn.RData")

# Raw counts
tb <- tibble(avg = rowMeans(counts(dds)), var = apply(counts(dds), 1, sd)^2)

ggplot(tb, aes(avg, var)) + geom_point() +
  geom_smooth() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format())) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format())) + 
  annotation_logticks() + labs(x = "Average", y = "Variance", 
                               title = "hESC project raw gene counts") +
  theme_bw() + theme(panel.grid.minor = element_blank())

# Normalised counts
sf <- log2(counts(dds, normalize = TRUE) + 1)
sf[, dds$rna_seq == "single_cell"] <- DrImpute::DrImpute(sf[, 
                                                              dds$rna_seq == 
                                                                "single_cell"])

tb <- tibble(avg = rowMeans(sf), var = apply(sf, 1, sd)^2)

ggplot(tb, aes(avg, var)) + geom_point() +
  geom_smooth() +
  labs(x = "Average", y = "Variance", 
       title = "hESC project normalised counts") +
  theme_bw() + theme(panel.grid.minor = element_blank())
