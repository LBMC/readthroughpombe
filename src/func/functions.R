#######################
##### VolcanoPlot #####
#######################
VolcanoPlot <- function(res_deseq, pval.thresh, where, save, title) {
  pdf(paste(where, save, sep = ""))
  plot(res_deseq$log2FoldChange, -log10(res_deseq$padj), xlab = "shrinked log2FC", ylab = "-log10(pval)", main = title) 
  nbup <- length(which(res_deseq$padj<pval.thresh & res_deseq$log2FoldChange>0))
  nbdown <- length(which(res_deseq$padj<pval.thresh & res_deseq$log2FoldChange<0))
  
  abline(v = 0, col = "gray", lty = 2)
  abline(h = -log10(pval.thresh), col = "red", lty = 2)     
  legend("topright", legend = paste("pval threshold=", pval.thresh, sep = ""), col = "red", lty = 2, cex =0.75, bty = "n")
  mtext(paste("nb up=", nbup, " - nb down=", nbdown, " on ", dim(res_deseq)[1], " genes", sep = ""), 3)
  dev.off()
}

####################
##### PvalHist #####
####################
PvalHist <- function(res_deseq, pval.thresh, where, save, title) {
  pdf(paste(where, save, sep = ""))
  hist(res_deseq$padj, breaks = 50, xlab = "pval", main = title )  
  nbupdown <- length(which(res_deseq$padj<pval.thresh))
  mtext(paste("nb diff=", nbupdown, " on ", dim(res_deseq)[1], " genes", sep = ""), 3)
  abline(v = pval.thresh, col = "red", lty = 2)     
  legend("topright", legend = paste("pval threshold=", pval.thresh, sep = ""), col = "red", lty = 2, cex =0.75, bty = "n")
  dev.off()
}

############################
##### MAplotBetween2Cn #####
############################
MAplotBetween2Cn <- function(res_deseq, where, save, title, alpha, hlim) {
  pdf(paste(where, save, sep = ""))
  plotMA(res_deseq, ylim=c(-3,3), main = title, alpha = alpha, MLE = F)
  abline(h = hlim, col = "dodgerblue", lwd = 2, lty = 2)
  dev.off()
}

#####################
##### HeatmapDE #####
#####################
HeatmapDE <- function(res_deseq, rld_count, pval.thresh, where, save, title) {
  ind.ord.pval <- order(res_deseq$padj)
  ind.signif <- which(res_deseq$padj[ind.ord.pval] <pval.thresh & is.na(res_deseq$padj[ind.ord.pval]) == F)
  de <- ind.signif[1:min(c(100, length(ind.signif)))]
  de.to.clust <-t(assay(rld_count)[rownames(res_deseq)[ind.ord.pval][de], ])
  pdf(paste(where, save, sep = ""))
  heatmap.2(de.to.clust, xlab = "Genes", ylab = "Samples", margin = c(8, 8), cexRow = 0.9, cex.main = 0.75, main = paste(title, " (rlog)", sep = ""), trace = "none",symm = F)
  dev.off()
}

##############################
##### DistBetweenSamples #####
##############################
DistBetweenSamples <- function(rld_count, where, save) {
  save.pref <- strsplit(save, ".", fixed = T)[[1]][1]
  pca <- plotPCA(rld_count, intgroup="condition")	
  pdf(paste(where, save.pref, "2.pdf", sep = ""))
  print(pca)
  dev.off()
}

########################
##### PlotCountDen #####
########################
PlotCountDen <- function(count_df, where, save) {
  den <- ggplot(count_df, aes(log2(value+1), fill = samples, linetype=replicate)) + geom_density(alpha = 0.5) + xlab("") + ylab(expression("density"(log[2](count + 1))))
  pdf(paste(where, save, sep = ""))
  print(den)
  dev.off()  
}

#############################
##### ComputeDistKSTest #####
#############################
ComputeDistKSTest <- function(dist1 = dist.diff, dist2 = dist.all) {
  cdf1 <- ecdf(dist1);  cdf2 <- ecdf(dist2)
  minMax <- seq(min(dist1, dist2), max(dist1, dist2), length.out=length(dist1)) 
  x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )] 
  return(x0)
}
