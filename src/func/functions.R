#######################  VolcanoPlot #####
VolcanoPlot <- function(res_deseq, pval.thresh, where, save, title) {
  pdf(paste(where, save, sep = ""))
  plot(res_deseq$log2FoldChange, -log10(res_deseq$padj), xlab = "shrinked log2FC", 
    ylab = "-log10(pval)", main = title)
  nbup <- length(which(res_deseq$padj < pval.thresh & res_deseq$log2FoldChange > 
    0))
  nbdown <- length(which(res_deseq$padj < pval.thresh & res_deseq$log2FoldChange < 
    0))
  
  abline(v = 0, col = "gray", lty = 2)
  abline(h = -log10(pval.thresh), col = "red", lty = 2)
  legend("topright", legend = paste("pval threshold=", pval.thresh, sep = ""), 
    col = "red", lty = 2, cex = 0.75, bty = "n")
  mtext(paste("nb up=", nbup, " - nb down=", nbdown, " on ", dim(res_deseq)[1], 
    " genes", sep = ""), 3)
  dev.off()
}

####################  PvalHist #####
PvalHist <- function(res_deseq, pval.thresh, where, save, title) {
  pdf(paste(where, save, sep = ""))
  hist(res_deseq$padj, breaks = 50, xlab = "pval", main = title)
  nbupdown <- length(which(res_deseq$padj < pval.thresh))
  mtext(paste("nb diff=", nbupdown, " on ", dim(res_deseq)[1], " genes", 
    sep = ""), 3)
  abline(v = pval.thresh, col = "red", lty = 2)
  legend("topright", legend = paste("pval threshold=", pval.thresh, sep = ""), 
    col = "red", lty = 2, cex = 0.75, bty = "n")
  dev.off()
}

############################ MAplotBetween2Cn #####
MAplotBetween2Cn <- function(res_deseq, where, save, title, pval.thresh, 
  hlim) {
  res.plot <- res_deseq[, c("baseMean", "log2FoldChange")]
  res.plot$signif <- res_deseq$padj < pval.thresh
  pdf(paste(where, save, sep = ""))
  plotMA(res.plot, ylim = c(-3, 3), main = title)
  abline(h = hlim, col = "dodgerblue", lwd = 2, lty = 2)
  dev.off()
}

##################### HeatmapDE #####
HeatmapDE <- function(res_deseq, rld_count, pval.thresh, where, save, title) {
  ind.ord.pval <- order(res_deseq$padj)
  ind.signif <- which(res_deseq$padj[ind.ord.pval] < pval.thresh & is.na(res_deseq$padj[ind.ord.pval]) == 
    F)
  de <- ind.signif[1:min(c(100, length(ind.signif)))]
  de.to.clust <- t(assay(rld_count)[rownames(res_deseq)[ind.ord.pval][de], 
    ])
  pdf(paste(where, save, sep = ""))
  heatmap.2(de.to.clust, xlab = "Genes", ylab = "Samples", margin = c(8, 
    8), cexRow = 0.9, cex.main = 0.75, main = paste(title, " (rlog)", 
    sep = ""), trace = "none", symm = F)
  dev.off()
}

############################## DistBetweenSamples #####
DistBetweenSamples <- function(rld_count, where, save) {
  save.pref <- strsplit(save, ".", fixed = T)[[1]][1]
  pca <- plotPCA(rld_count, intgroup = "condition")
  pdf(paste(where, save.pref, "2.pdf", sep = ""))
  print(pca)
  dev.off()
}

######################## PlotCountDen #####
PlotCountDen <- function(count_df, where, save) {
  den <- ggplot(count_df, aes(log2(value + 1), fill = samples, linetype = replicate)) + 
    geom_density(alpha = 0.5) + xlab("") + ylab(expression(density(log[2](count + 
    1))))
  pdf(paste(where, save, sep = ""))
  print(den)
  dev.off()
}

############################# ComputeDistKSTest #####
ComputeDistKSTest <- function(dist1 = dist.diff, dist2 = dist.all) {
  cdf1 <- ecdf(dist1)
  cdf2 <- ecdf(dist2)
  minMax <- seq(min(dist1, dist2), max(dist1, dist2), length.out = length(dist1))
  x0 <- minMax[which(abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - 
    cdf2(minMax))))]
  return(x0)
}

###################### ViolinPlot #####
ViolinPlot <- function(input.df.with.count, ylabel, save.pdf, do.col = F) {
  if (do.col) {
    input.df.with.count$genes <- as.factor(input.df.with.count$genes)
    n <- length(levels(input.df.with.count$genes))
    p <- ggplot(input.df.with.count, aes(factor(samples), log2(value + 
      1), alpha = 0.3)) + geom_violin(draw_quantiles = c(0.25, 0.5, 
      0.75)) + theme(legend.position = "none") + geom_jitter(aes(colour = genes), 
      height = 0, width = 0.2) + facet_grid(~chromosome) + ylab(ylabel) + 
      xlab("condition") + scale_color_manual(values = rainbow(n))
  } else {
    p <- ggplot(input.df.with.count, aes(factor(samples), log2(value + 
      1), alpha = 0.3)) + geom_violin(draw_quantiles = c(0.25, 0.5, 
      0.75)) + theme(legend.position = "none") + geom_jitter(height = 0, 
      width = 0.2) + facet_grid(~chromosome) + ylab(ylabel) + xlab("condition")
  }
  pdf(save.pdf)
  print(p)
  dev.off()
}

################################### ComputeWindowsOverChrom #####
ComputeWindowsOverChrom <- function(vect.pos, size.chr) {
  tab <- table(vect.pos)
  seq.along.chr <- 1:size.chr
  missing.pos <- setdiff(seq.along.chr, as.numeric(names(tab)))
  tab.missing <- rep(0, length(missing.pos))
  names(tab.missing) <- missing.pos
  tab.all <- c(tab, tab.missing)
  tab.all <- tab.all[order(as.numeric(names(tab.all)))]
  return(tab.all)
}
