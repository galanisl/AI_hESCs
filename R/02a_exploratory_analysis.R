library(DESeq2)
library(readr)
library(dplyr)
library(cowplot)
library(ggplot2)
source("utility_functions.R")

# Data preparation --------------------------------------------------------

# Load the count data
load("results/hesc_all.RData")

# Read in the sample aliases and arrange to match the count matrix column order
stype <- read_tsv("all_samples.tsv")
stype <- stype %>% 
  mutate(sample_type = factor(paste0(cell_line, " ", ifelse(medium == "none", 
                                                            "", medium), " (", 
                                     stringr::str_to_title(batch), " et al.)"),
                              levels = c("Epi  (Blakeley et al.)", 
                                         "Epi  (Petropoulos et al.)",
                                         "Epi  (Yan et al.)",
                                         "PE  (Blakeley et al.)", 
                                         "PE  (Petropoulos et al.)",
                                         "TE  (Blakeley et al.)",
                                         "CH1 AI (Wamaitha et al.)",
                                         "CH2 AI (Wamaitha et al.)",
                                         "CH3 AI (Wamaitha et al.)",
                                         "H1 AI (Wamaitha et al.)",
                                         "H9 AI (Wamaitha et al.)",
                                         "H9 mTeSR1_matrigel (Wamaitha et al.)",
                                         "H1 mTeSR1_laminin (Wamaitha et al.)",
                                         "H9 mTeSR1_laminin (Wamaitha et al.)",
                                         "Derived_p0 KSR/FBS+FGF (Yan et al.)",
                                         "Derived_p10 KSR/FBS+FGF+L (Yan et al.)",
                                         "WIBR3 KSR/FBS+FGF (Ji et al.)",
                                         "H9 KSR/FGF (Takashima et al.)",
                                         "Shef6 KSR/FGF (Guo et al.)",
                                         "WIBR3 5iLA (Ji et al.)",
                                         "H9 t2iL+Go (Guo et al.)",
                                         "H9 t2iL+Go (Takashima et al.)",
                                         "HNES1 t2iL+Go (Guo et al.)",
                                         "Shef6_p18 t2iL+Go (Guo et al.)",
                                         "Shef6_p26 t2iL+Go (Guo et al.)",
                                         "H9 t2iL+Go+Y (Guo et al.)",
                                         "HNES1 t2iL+Go FF (Guo et al.)",
                                         "Shef6 t2iL+Go FF (Guo et al.)"), 
                              ordered = TRUE))
stype <- stype %>% 
  mutate(condition = parse_factor(condition, NULL), 
         exp_group = parse_factor(exp_group, NULL),
         rna_seq = parse_factor(rna_seq, NULL),
         cell_line = parse_factor(cell_line, NULL),
         medium = parse_factor(medium, NULL))

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

# Putting all the data together
hesc$abundance <- hesc$abundance[, stype$sample_name]
hesc$counts <- hesc$counts[, stype$sample_name]
hesc$length <- hesc$length[, stype$sample_name]
dds <- DESeqDataSetFromTximport(txi = hesc, colData = stype, 
                                design = ~ exp_group)

# Removal of mitochondrial-, ribosomal- and pseudo-genes
mito <- grep("^MT-", rownames(dds))
dds <- dds[-mito, ]

ribo <- read_tsv("data/ribosomal.txt")
dds <- dds[!(rownames(dds) %in% ribo$`Approved Symbol`), ]

load("data/t2g_complete.RData")

pseudo <- t2g[grep("pseudogene", t2g$gene_biotype), ]
dds <- dds[!(rownames(dds) %in% pseudo$symbol), ]

# Removal of no-show genes
dds <- dds[!near(rowSums(counts(dds)), 0), ]

# Estimate size factors
dds <- estimateSizeFactors(dds)

# Normalise and impute the single cells
sf <- log2(counts(dds, normalize = TRUE) + 1)
sf[, stype$rna_seq == "single_cell"] <- DrImpute::DrImpute(sf[, 
                                                              stype$rna_seq == 
                                                                "single_cell"])

# Removal of invariant genes
meanDM <- mean(sf)
nSD <- apply(sf, 1, function(x) sd(x)/meanDM)
if(sum(nSD < 1e-3) > 0){
  sf <- sf[-which(nSD < 1e-3), ]
  dds <- dds[-which(nSD < 1e-3), ]
}

# PCA before batch effect removal -----------------------------------------

# p-value for variable genes
pval <- 0.01

pca_wbe <- perform_pca(sf, pval = pval, top_pc = 5)
colnames(pca_wbe$ind$coord) <- paste0("PC", 1:5)

# Pairs plot coloured by RNA-seq type
pairs(pca_wbe$ind$coord, col = factor(stype$rna_seq), 
      main = paste0("Before batch effect removal (n = ", 
                    nrow(pca_wbe$var$coord), ", adjusted P = ", pval, ")"), 
      pch = 16, oma = c(3, 3, 7, 15), cex = 0.7)
legend("right", fill = unique(factor(stype$rna_seq)), 
       legend = unique(stype$rna_seq))

# Batch effect removal ----------------------------------------------------

# Removal of batch effects
sf <- limma::removeBatchEffect(sf, batch = stype$rna_seq)

pca_nbe <- perform_pca(sf, pval = pval, top_pc = 5)
colnames(pca_nbe$ind$coord) <- paste0("PC", 1:5)

# Pairs plot coloured by RNA-seq type
pairs(pca_nbe$ind$coord, col = factor(stype$rna_seq), 
      main = paste0("After batch effect removal (n = ", 
                    nrow(pca_nbe$var$coord), ", adjusted P = ", pval, ")"), 
      pch = 16, oma = c(3, 3, 7, 15), cex = 0.7)
legend("right", fill = unique(factor(stype$rna_seq)), 
       legend = unique(stype$rna_seq))

p <- list()
p[[1]] <- plot_pca(pca_nbe, stype$sample_type, 1, 2, pval, col_pal, shape_pal,
                   "After batch effect removal") + theme(legend.position = "none")
p[[2]] <- NULL
p[[3]] <- plot_pca(pca_nbe, stype$sample_type, 1, 3, pval, col_pal, shape_pal,
                   "After batch effect removal") + theme(legend.position = "none")
p[[4]] <- plot_pca(pca_nbe, stype$sample_type, 2, 3, pval, col_pal, shape_pal,
                   "After batch effect removal") + theme(legend.position = "none")
plot_grid(plotlist = p, nrow = 2, ncol = 2, align = "hv")

save(pca_wbe, pca_nbe, col_pal, shape_pal, pval, file = "results/pca_sfn.RData")
save(dds, sf, file = "results/raw_norm_DESeq_sfn.RData")

# 3D PCA ------------------------------------------------------------------
library(plotly)

res <- tibble(x = pca_nbe$ind$coord[, 1], y = pca_nbe$ind$coord[, 2], 
              z = pca_nbe$ind$coord[, 3],
              cond = stype$sample_type)

pc1_var <- round(pca_nbe$eig[1, "percentage of variance"],1)
pc2_var <- round(pca_nbe$eig[2, "percentage of variance"],1)
pc3_var <- round(pca_nbe$eig[3, "percentage of variance"],1)

plot_ly(res, x = ~x, y = ~y, z = ~z, color = ~cond, marker = list(size = 5),
        colors = col_pal) %>%
  add_markers() %>%
  layout(title = paste0("After batch effect removal (n = ", 
                        nrow(pca_nbe$var$coord), ", adjusted P = ", pval, ")"),
         scene = list(xaxis = list(title = paste0("PC1 (", pc1_var, "%)")),
                      yaxis = list(title = paste0("PC2 (", pc2_var, "%)")),
                      zaxis = list(title = paste0("PC3 (", pc3_var, "%)"))))

# tSNE --------------------------------------------------------------------
library(Rtsne)

# Run 100 t-SNEs and keep the one with the lowest Kullback-Leibler divergence
# (recommended by Laurens van der Maaten in https://lvdmaaten.github.io/tsne/)
fit <- scran::trendVar(sf)
decomp <- scran::decomposeVar(sf, fit)
idx <- which(decomp$FDR < pval)
ntop <- length(idx)
mat <- t(sf[idx,])

kl_div <- Inf
for(i in 1:100){
  ts <- Rtsne(mat, dims = 3, perplexity = 30)
  if(min(ts$itercosts) < kl_div){
    kl_div <- min(ts$itercosts)
    tsf <- ts
  }
}

p <- list()
p[[1]] <- plot_tsne(tsf, stype$sample_type, d1 = 1, d2 = 2, pval = pval, ntop, 
                    col_pal, shape_pal, "After batch effect removal") + 
  theme(legend.position = "none")
p[[2]] <- NULL
p[[3]] <- plot_tsne(tsf, stype$sample_type, d1 = 1, d2 = 3, pval = pval, ntop, 
                    col_pal, shape_pal, "After batch effect removal") + 
  theme(legend.position = "none")
p[[4]] <- plot_tsne(tsf, stype$sample_type, d1 = 2, d2 = 3, pval = pval, ntop, 
                    col_pal, shape_pal, "After batch effect removal") + 
  theme(legend.position = "none")
plot_grid(plotlist = p, nrow = 2, ncol = 2, align = "hv")

res_ts <- tibble(x = tsf$Y[, 1], y = tsf$Y[, 2], z = tsf$Y[, 3],
                 cond = stype$sample_type)

plot_ly(res_ts, x = ~x, y = ~y, z = ~z, color = ~cond, marker = list(size = 5),
        colors = col_pal) %>%
  add_markers() %>%
  layout(title = paste0("After batch effect removal (n = ", ntop, 
                        ", adjusted P = ", pval, ", 100 runs)"),
         scene = list(xaxis = list(title = "t-SNE1"),
                      yaxis = list(title = "t-SNE2"),
                      zaxis = list(title = "t-SNE3")))

save(tsf, res_ts, col_pal, shape_pal, pval, ntop, file = "results/tsne_sfn.RData")
var_genes <- colnames(mat)
save(var_genes, file = "results/var_genes.RData")


# UMAP --------------------------------------------------------------------

usf <- uwot::umap(mat, pca = 50, n_components = 3)

p <- list()
p[[1]] <- plot_umap(usf, stype$sample_type, d1 = 1, d2 = 2, pval, ntop, 
                    col_pal, shape_pal) + 
  theme(legend.position = "none")
p[[2]] <- NULL
p[[3]] <- plot_umap(usf, stype$sample_type, d1 = 1, d2 = 3, pval, ntop,  
                    col_pal, shape_pal) + 
  theme(legend.position = "none")
p[[4]] <- plot_umap(usf, stype$sample_type, d1 = 2, d2 = 3, pval, ntop,  
                    col_pal, shape_pal) + 
  theme(legend.position = "none")
plot_grid(plotlist = p, nrow = 2, ncol = 2, align = "hv")

save(usf, col_pal, shape_pal, ntop, pval, file = "results/umap_sfn.RData")