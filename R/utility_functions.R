
# Utility functions -------------------------------------------------------

# Perform PCA of X using the most variable genes based on pval
perform_pca <- function(X, pval = 0.001, top_pc = 5){
  # Detecting variable genes (corrected p-val < 0.001)
  fit <- scran::trendVar(X)
  decomp <- scran::decomposeVar(X, fit)
  idx <- which(decomp$FDR < pval)
  ntop <- length(idx)
  
  mat <- t(X[idx, ])
  
  pca <- FactoMineR::PCA(mat, scale.unit = FALSE, ncp = top_pc, graph = FALSE)
  return(pca)
}

# Plot 2D PCA for the given pca object, the given components pc1, pc2, and the
# given colour palette and title. Points are coloured by condition cond.
plot_pca <- function(pca, cond, pc1 = 1, pc2 = 2, pval = 0.001, col_pal, shape_pal, ptitle){
  res <- tibble(x = pca$ind$coord[, pc1], y = pca$ind$coord[, pc2], 
                cond = cond)
  
  pc1_var <- round(pca$eig[pc1, "percentage of variance"],1)
  pc2_var <- round(pca$eig[pc2, "percentage of variance"],1)
  
  p_pca <- ggplot(res, aes(x, y, colour = cond, shape = cond)) + 
    geom_point(size = 2.5) +
    scale_colour_manual(values = col_pal) + 
    scale_shape_manual(values = shape_pal) + 
    labs(x = paste0("PC", pc1, " (", pc1_var, "%)"), 
         y = paste0("PC", pc2, " (", pc2_var, "%)"),
         title = paste0(ptitle, " (n = ", nrow(pca$var$coord),
                        ", adjusted P = ", pval, ")"), 
         colour = "", shape = "") + 
    theme_bw(base_size = 15)
  return(p_pca)
}

# Plot 2D t-SNE for the given t-SNE object, the given dimensions d1, d2, and the
# given colour palette and title. Points are coloured by condition cond.
plot_tsne <- function(tsf, cond, d1 = 1, d2 = 2, pval = 0.001, ntop, col_pal, shape_pal, ptitle){

  res_ts <- tibble(x = tsf$Y[, d1], y = tsf$Y[, d2], cond = cond)
  
  p_tsne <- ggplot(res_ts, aes(x, y, colour = cond, shape = cond)) + 
    geom_point(size = 2.5) +
    scale_colour_manual(values = col_pal) + 
    scale_shape_manual(values = shape_pal) + 
    labs(x = paste0("t-SNE", d1), 
         y = paste0("t-SNE", d2), 
         title = paste0("After batch effect removal (n = ", ntop, 
                        ", adjusted P = ", pval, ", 100 runs)"), 
         colour = "", shape = "") + 
    theme_bw(base_size = 15)
  return(p_tsne)
}

# Plot 2D UMAP for the given UMAP object, the given dimensions d1, d2, and the
# given colour palette and title. Points are coloured by condition cond.
plot_umap <- function(usf, cond, d1 = 1, d2 = 2, ntop, col_pal, shape_pal, ptitle){
  
  res_us <- tibble(x = usf[, d1], y = usf[, d2], cond = cond)
  
  p_umap <- ggplot(res_us, aes(x, y, colour = cond, shape = cond)) + 
    geom_point(size = 2.5) +
    scale_colour_manual(values = col_pal) + 
    scale_shape_manual(values = shape_pal) + 
    labs(x = paste0("UMAP", d1), 
         y = paste0("UMAP", d2), 
         title = paste0("After batch effect removal (n = ", ntop, ")"), 
         colour = "", shape = "") + 
    theme_bw(base_size = 15)
  return(p_umap)
}

# Check how the most informative genes map onto PC1 and PC2
plot_bi <- function(pca, cond, pc1 = 1, pc2 = 2, top_genes = 100,
                    pval = 0.001, col_pal, ptitle){
  
  #u <- sweep(pca$ind$coord[, c(pc1, pc2)], 2, pca$svd$vs[c(1, 2)], FUN='*')
  u <- pca$ind$coord[, c(pc1, pc2)]
  
  res <- tibble(x = u[, 1], y = u[, 2], 
                cond = cond)
  # cos2 <- pca$var$cos2[, pc1] + pca$var$cos2[, pc2]
  cos2 <- pca$var$contrib[, pc1]
  cos2 <- sort(cos2, decreasing = TRUE, index.return = TRUE)
  
  # Correaltion circle
  r <- sqrt(qchisq(0.69, df = 2)) * prod(colMeans(u^2))^(1/4)
  
  # Scale directions
  v.scale <- rowSums(pca$var$coord^2)
  v <- r * pca$var$coord / sqrt(max(v.scale))
  
  genes <- tibble(x = v[cos2$ix[1:top_genes], pc1],
                  y = v[cos2$ix[1:top_genes], pc2],
                  varname = rownames(pca$var$coord[cos2$ix[1:top_genes], ]),
                  angle = (180/pi) * atan(y / x),
                  hjust = (1 - 1.1 * sign(x)) / 2)
  
  pc1_var <- round(pca$eig[pc1, "percentage of variance"],1)
  pc2_var <- round(pca$eig[pc2, "percentage of variance"],1)
  
  p_bi <- ggplot(res, aes(x, y, colour = cond)) + 
    geom_point(size = 2.5) +
    scale_colour_manual(values = col_pal) + 
    labs(x = paste0("PC", pc1, " (", pc1_var, "%)"), 
         y = paste0("PC", pc2, " (", pc2_var, "%)"),
         title = paste0(ptitle, " (n = ", nrow(pca$var$coord),
                        ", adjusted P = ", pval, ")"), 
         colour = "") + 
    theme_bw(base_size = 15)
  p_bi <- p_bi + 
    geom_segment(data = genes, aes(x = 0, y = 0, xend = x, yend = y),
                 arrow = arrow(length = unit(1/2, "picas")), 
                 color = "darkred") +
    geom_text(data = genes, aes(label = varname, x = x, y = y, 
                                angle = angle, hjust = hjust), 
              color = "darkred", size = 3)
  return(p_bi)
  
}

#' This is the function for selection of overdispersed genes adapted from:
#' https://github.com/10XGenomics/single-cell-3prime-paper/
#'
#' @param m  a (protentially sparse) gene x cells count matrix
#' @return a vector of normalized (robust Z-scores) dispersion values, 
#' one per gene.
select_variable_genes <- function(m) {
  df <- data.frame(mean = rowMeans(m + 1/ncol(m)), 
                   cv = apply(m,1,sd) / rowMeans(m + 1/ncol(m)), 
                   var = apply(m,1,var))
  df$dispersion <- with(df, var/mean)
  df$mean_bin <- with(df, 
                      cut(mean, breaks = c(-Inf, 
                                           unique(quantile(mean, 
                                                           seq(0.1,1,0.05), 
                                                           na.rm = TRUE)), 
                                           Inf)))
  var_by_bin <- data.frame(
    mean_bin = factor(levels(df$mean_bin), levels = levels(df$mean_bin)),
    bin_median = as.numeric(tapply(df$dispersion, df$mean_bin, stats::median)),
    bin_mad = as.numeric(tapply(df$dispersion, 
                                df$mean_bin, 
                                stats::mad)))[table(df$mean_bin) > 0,]
  df$bin_disp_median <- var_by_bin$bin_median[match(df$mean_bin, 
                                                    var_by_bin$mean_bin)]
  df$bin_disp_mad <- var_by_bin$bin_mad[match(df$mean_bin, 
                                              var_by_bin$mean_bin)]
  df$dispersion_norm <- with(
    df, (dispersion - bin_disp_median)/(bin_disp_mad + 0.01)
  )
  return(df$dispersion_norm)
}

# Transform DESeqResults into a tibble
DESeqRes2tibble <- function(dres){
  dres$gene <- rownames(dres)
  
  tb <- tibble(gene = dres$gene, log2FoldChange = dres$log2FoldChange, 
               pvalue = dres$pvalue, padj = dres$padj, baseMean = dres$baseMean, 
               lfcSE = dres$lfcSE, stat = dres$stat) %>% 
    arrange(desc(log2FoldChange), padj)

  return(tb)
}

# Volcano plot from DESeqResults converted to tibble
plot_volcano <- function(dtb, pval = 0.05, log2fc = 1, ptitle){
  v <- ggplot(dtb, aes(log2FoldChange, -log10(padj), label = gene)) +
    geom_point(shape = 1, alpha = 0.2) + 
    geom_text_repel(data = filter(dtb, (padj <= pval & abs(log2FoldChange) >= log2fc)),
                    colour = "blue") +
    geom_hline(yintercept = -log10(pval), linetype = 2, colour = "red") +
    geom_vline(xintercept = -log2fc, linetype = 2, colour = "red") +
    geom_vline(xintercept = log2fc, linetype = 2, colour = "red") +
    labs(x = expression(paste(log[2], "-fold-change")), 
         y = expression(paste(-log[10], "(p-value)")),
         title = ptitle) +
    theme_bw(base_size = 15)
  return(v)
}

# Given a DESeqResults converted to tibble and log2FC/p-value thresholds, it
# performs a functional enrichment analysis of the up- and down-regulated genes.
# The 'benj' parameter indicates whether adjusted or unadjusted p-values are 
# plotted
funenrich_analysis <- function(deseq_tb, log2FC, pval, benj = TRUE){
  up <- deseq_tb %>% 
    filter(log2FoldChange >= log2FC & padj <= pval)
  dw <- deseq_tb %>% 
    filter(log2FoldChange <= -log2FC & padj <= pval)
  
  enr_up <- fun_enrich(up$gene, deseq_tb$gene, "SYMBOL", benjamini = TRUE)
  penr_up <- plot_fun_enrich(enr_up, benjamini = benj, char_per_line = 60) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "red") +
    theme_bw(base_size = 12) + theme(legend.title = element_blank(), 
                                     legend.background = element_blank(), 
                                     legend.position="top")
  
  enr_dw <- fun_enrich(dw$gene, deseq_tb$gene, "SYMBOL", benjamini = TRUE)
  penr_dw <- plot_fun_enrich(enr_dw, benjamini = benj, char_per_line = 60) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "red") +
    theme_bw(base_size = 12) + theme(legend.title = element_blank(), 
                                     legend.background = element_blank(), 
                                     legend.position="top")
  
  return(list(up = up, dw = dw, enr_up = enr_up, enr_dw = enr_dw, 
              p = plot_grid(penr_up, penr_dw, nrow = 1, ncol = 2)))
}
