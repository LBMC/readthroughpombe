require(sleuth)
require(rhdf5)
require(stringr)
require(dplyr)
require(data.table)
require(zoo)
require(MASS)
require(VennDiagram)
source("src/func/functions.R")

wt <- paste("wt", 1:3, sep = "")
cut14 <- paste("cut14", 1:3, sep = "")
dble <- paste("cut14_cdc15", 1:3, sep = "")
cdc15 <- paste("cdc15", 1:3, sep = "")
rrp6D <- paste("rrp6D", 1:3, sep = "")
alpha_thresh <- 0.01
do.perm.dist <- F

##### Read annotation file
##### Load annotations in bed 
bed <- fread("data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr.bed", 
  h = F, sep = "\t", stringsAsFactors = F)
bed.transcripts <- as.data.frame(bed[which(bed$V8 == "transcript" & bed$V1 %in% c("I", "II", "III")), ])
bed.transcripts$V2 <- bed.transcripts$V2+1 #only use for plots
bed.transcripts$deb <- apply(bed.transcripts, 1, function(x) if (x["V6"] == "+") {
  x["V2"]
} else {
  x["V3"]
})
bed.transcripts$deb <- as.numeric(bed.transcripts$deb)

##### tRNA
tRNA.index <- sapply(bed.transcripts$V10, function(x) str_detect(x, "tRNA"))
bed.tRNA <- bed.transcripts[which(tRNA.index == T), ]
bed.transcripts <- bed.transcripts[-which(tRNA.index == T), ]

##### Sizes of chromosomes
sizes <- matrix(c(5579133, 4539804, 2452883, 19431,20218), nrow = 1, dimnames = list(NULL, c("I", 
  "II", "III", "MT", "MTR")))

##### Names of features for analysis without readthrough
names.all <- bed.transcripts$V4
##### Names of features for analysis with readthrough
bed_rrp6D_based <- fread("results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_rrp6D_wt_based_cut14_forward_reverse.bed", sep = "\t", h = F)
bed_rrp6D_based$V2 <- bed_rrp6D_based$V2+1 #only use for plots
bed_rrp6D_based$deb <- apply(bed_rrp6D_based, 1, function(x) if (x["V6"] == "+") {
  x["V2"]
} else {
  x["V3"]
})
bed_rrp6D_based$deb <- as.numeric(bed_rrp6D_based$deb)
bed_rrp6D_based <- bed_rrp6D_based[which(bed_rrp6D_based$V1 %in% c("I", "II", "III")), ]
names_bed_rrp6D_based <- unique(c(bed_rrp6D_based$V4, names.all))
names_bed_rrp6D_based_rt <- fread("results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_cut14_forward_reverse.bed", sep = "\t", 
  h = F)
names_bed_rrp6D_based_rt <- names_bed_rrp6D_based_rt[which(names_bed_rrp6D_based_rt$V1 %in% c("I", "II", "III")), ]$V4

bed_cut14_based <- fread("results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.bed", sep = "\t", h = F)
bed_cut14_based$V2 <- bed_cut14_based$V2+1 #only use for plots
bed_cut14_based$deb <- apply(bed_cut14_based, 1, function(x) if (x["V6"] == "+") {
  x["V2"]
} else {
  x["V3"]
})
bed_cut14_based$deb <- as.numeric(bed_cut14_based$deb)
bed_cut14_based <- bed_cut14_based[which(bed_cut14_based$V1 %in% c("I", "II", "III")), ]
names_bed_cut14_based <- unique(c(bed_cut14_based$V4, names.all))
names_bed_cut14_based_rt <- fread("results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_rrp6D_forward_reverse.bed", sep = "\t", 
  h = F)
names_bed_cut14_based_rt <- names_bed_cut14_based_rt[which(names_bed_cut14_based_rt$V1 %in% c("I", "II", "III")), ]$V4

##### Describe data to prepare sleuth input design
s2c.all.analysis <- list(data.frame(sample = c(wt, cut14), condition = c(rep("wt", 
  3), rep("cut14", 3)), path = c(paste("results/kallisto_quantif/wt_all_genes_R", 
  1:3, sep = ""), paste("results/kallisto_quantif/cut14_all_genes_R", 1:3, sep = "")), 
  stringsAsFactors = F), 
  data.frame(sample = c(wt, dble), condition = c(rep("wt", 
  3), rep("cut14_cdc15", 3)), path = c(paste("results/kallisto_quantif/wt_all_genes_R", 
  1:3, sep = ""), paste("results/kallisto_quantif/cut14_cdc15_all_genes_R", 1:3, 
  sep = "")), stringsAsFactors = F), 
  data.frame(sample = c(wt, rrp6D), condition = c(rep("wt", 
  3), rep("rrp6D", 3)), path = c(paste("results/kallisto_quantif/wt_all_genes_R", 
  1:3, sep = ""), paste("results/kallisto_quantif/rrp6D_all_genes_R", 1:3, sep = "")), 
  stringsAsFactors = F),
  data.frame(sample = c(wt, cdc15), condition = c(rep("wt", 
  3), rep("cdc15", 3)), path = c(paste("results/kallisto_quantif/wt_all_genes_R", 
  1:3, sep = ""), paste("results/kallisto_quantif/cdc15_all_genes_R", 1:3, sep = "")), 
  stringsAsFactors = F),
 
  data.frame(sample = c(wt, cut14), condition = c(rep("wt", 
  3), rep("cut14", 3)), path = c(paste("results/kallisto_quantif/wt_cut14wt_based_+rrp6D_readthrough_with_other_genes_R", 
  1:3, sep = ""), paste("results/kallisto_quantif/cut14_cut14wt_based_+rrp6D_readthrough_with_other_genes_R", 1:3, sep = "")), 
  stringsAsFactors = F), 
  data.frame(sample = c(wt, dble), condition = c(rep("wt", 
  3), rep("cut14_cdc15", 3)), path = c(paste("results/kallisto_quantif/wt_cut14wt_based_+rrp6D_readthrough_with_other_genes_R", 
  1:3, sep = ""), paste("results/kallisto_quantif/cut14_cdc15_cut14wt_based_+rrp6D_readthrough_with_other_genes_R", 1:3, 
  sep = "")), stringsAsFactors = F), 
  data.frame(sample = c(wt, rrp6D), condition = c(rep("wt", 
  3), rep("rrp6D", 3)), path = c(paste("results/kallisto_quantif/wt_rrp6Dwt_based_+cut14_readthrough_with_other_genes_R", 
  1:3, sep = ""), paste("results/kallisto_quantif/rrp6D_rrp6Dwt_based_+cut14_readthrough_with_other_genes_R", 1:3, sep = "")), 
  stringsAsFactors = F))

names(s2c.all.analysis) <- c("cut14_wt", "cut14_cdc15_wt", "rrp6D_wt",  "cdc15_wt", "cut14_wt_rt_only", "cut14_cdc15_wt_rt_only", "rrp6D_wt_rt_only")
bed.used <- list(names.all, names.all, names.all, names.all, names_bed_cut14_based_rt, names_bed_cut14_based_rt, names_bed_rrp6D_based_rt)
names(bed.used) <- c("cut14_wt", "cut14_wt", "rrp6D_wt", "cdc15_wt", "cut14_wt", "cut14_wt", "rrp6D_wt")
# Define ref and testing level for the strain
ref.level <- rep("wt", length(s2c.all.analysis))
comp.level <- c("cut14", "cut14_cdc15", "rrp6D", "cdc15", "cut14", "cut14_cdc15", "rrp6D")


##### Subset from kallisto
# For the purposes of analyses, need to fiter the quantification files produces by kallisto created new tsv and new subsetted h5 files: this thread was helpful https://github.com/pachterlab/sleuth/issues/131
fname <- "abundance.h5"
for (j in seq_along(s2c.all.analysis)){
  keep <- bed.used[[j]]
  s2c.tmp <- s2c.all.analysis[[j]]
  for (path in unique(s2c.tmp$path)) {
    suff <- ""
    if(str_detect(names(s2c.all.analysis)[j], "_rt_only")){
      suff <- "_rt"
    }
    system(paste("mkdir ", path, "_subset", suff, sep = ""))     

    all <- read_kallisto(path)
    file.remove(fname)
    name <- paste(path, "_subset", suff, "/", fname, sep = "")
    subset <- which(all$abundance$target_id %in% keep)  # all + rt
    rhdf5::h5createFile(fname)
    rhdf5::h5createGroup(fname, "aux")
    dims <- c(length(subset), 1)
    rhdf5::h5createDataset(fname, "aux/ids", dims = dims, storage.mode = "character", 
      size = 100, level = 7L)
    rhdf5::h5write(all$abundance$target_id[subset], fname, "aux/ids")
    new_boot <- lapply(all$bootstrap, function(x) x[subset, ])
    rhdf5::h5write(length(new_boot), fname, "aux/num_bootstrap")
    rhdf5::h5createGroup(fname, "bootstrap")
    for (i in seq_along(new_boot)) {
      bs <- new_boot[[i]]$est_counts
      bs_name <- paste0("bootstrap/bs", i - 1)
      rhdf5::h5write(bs, fname, bs_name)
    }
    rhdf5::h5write(all$abundance$eff_len[subset], fname, "aux/eff_lengths")
    rhdf5::h5write(all$abundance$len[subset], fname, "aux/lengths")
    rhdf5::h5write(all$abundance$est_counts[subset], fname, "est_counts")
    rhdf5::h5write(all$fld[subset], fname, "aux/fld")
    rhdf5::h5write(all$bias_normalized[subset], fname, "aux/bias_normalized")
    rhdf5::h5write(all$bias_observed[subset], fname, "aux/bias_observed")
    rhdf5::h5write(attributes(all)$num_processed, fname, "aux/num_processed")
    rhdf5::h5write(attributes(all)$index_version, fname, "aux/index_version")
    rhdf5::h5write(attributes(all)$kallisto_version, fname, "aux/kallisto_version")
    rhdf5::h5write(attributes(all)$start_time, fname, "aux/start_time")
    system(paste("mv ", fname, " ", name, sep = ""))
    H5close()
  }
}


##### sleuth analysis
for (i in seq_along(s2c.all.analysis)){
  save.path.dir <- paste("results/readthrough_analysis/sleuth_analysis/", sep = "")
  save.dir <- names(s2c.all.analysis)[i] 
  # Create results directory
  system(paste("mkdir ", save.path.dir, save.dir, sep = ""))
  s2c <- s2c.all.analysis[[i]]
  suff <- ""
  if(str_detect(names(s2c.all.analysis)[i], "_rt_only")){
    suff <- "_rt"
  }
  s2c$path <- paste(s2c$path, paste("_subset", suff, sep = ""), sep = "")
  s2c$condition <- as.factor(s2c$condition)
  s2c$condition = relevel(s2c$condition, ref= ref.level[i])
  code <- matrix(s2c$condition, dimnames = list(NULL, s2c$sample), nrow = 1)
  # Prepare sleuth input
  so <- sleuth_prep(s2c, extra_bootstrap_summary = T)
  # PCA plot 
  pca <- plot_pca(so, color_by = 'condition', text_labels = TRUE)
  pdf(paste(save.path.dir, save.dir, "/PCAplot_", save.dir, ".pdf", sep = ""))
  print(pca)
  dev.off()
  pca_axis <-   plot_pc_variance(so, use_filtered = T, pca_number = 3, scale = T)
  pdf(paste(save.path.dir, save.dir, "/PCAAxisplot_", save.dir, ".pdf", sep = ""))
  print(pca_axis)
  dev.off()

  if(class(so) != "try-error") {
    so <- sleuth_fit(so, ~condition, 'full')
    so <- sleuth_fit(so, ~1, 'reduced')

    #so <- sleuth_lrt(so, 'reduced', 'full')
    #sleuth_table_lrt <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

    ##### Are there differentially expressed genes?
    comp <- paste("condition", comp.level[i], sep = "")
    so <- sleuth_wt(so, comp, 'full') # change within loop
    sleuth_table_wt <- sleuth_results(so, comp, ref.level[i], show_all = T)
    sleuth_table_wt_filt <- sleuth_results(so, comp, ref.level[i], show_all = F)
    #show_all: if ‘TRUE’ will show all transcripts (not only the ones
    #  passing filters). The transcripts that do not pass filters
    #  will have ‘NA’ values in most columns.
    #pval: p-value of the chosen model
    #qval: false discovery rate adjusted p-value, using BH
    #b: 'beta' value (effect size). Technically a biased estimator of the fold change
    #se_b: standard error of the beta
    #mean_obs: mean of natural log counts of observations
    #var_obs: variance of observation
    #tech_var: technical variance of observation from the bootstraps
    #sigma_sq: raw estimator of the variance once the technical variance has been removed
    #smooth_sigma_sq: smooth regression fit for the shrinkage estimation
    #final_simga_sq: max(sigma_sq, smooth_sigma_sq); used forcovariance estimation of beta
    #sleuth_table_wt_filt$ratio <- exp(sleuth_table_wt_filt$b)
    #sleuth_table_wt_filt$proxy_log2_ratio <- log2(exp(sleuth_table_wt_filt$b))
    init <- dim( sleuth_results(so, comp, ref.level[i], show_all = T))[1]
    filt <- dim(sleuth_table_wt_filt)[1]
    sleuth_table_wt_filt$log2FC <-log2( exp(sleuth_table_wt_filt$b))
    sleuth_table_wt$log2FC <- log2( exp(sleuth_table_wt$b))
    write.table(sleuth_table_wt, file = paste(save.path.dir, save.dir, "/sleuth_shrinked_analysis_diff_expressed_", save.dir, ".txt", sep = ""), sep = "\t", col.names = T, row.names = T, quote = F)

    ##### VolcanoPlot 
    volcano <- plot_volcano(so, comp, test_type = "wt", which_model = "full", sig_level = alpha_thresh, point_alpha = 0.2, sig_color = "red", highlight = NULL)
    pdf(paste(save.path.dir, save.dir, "/volcano_", save.dir, ".pdf", sep = ""))
    print(volcano+ggtitle(paste("nb init = ", init, " - after sleuth filt = ", filt, sep = "")))
    dev.off()
		
    ##### Pval hist
    pdf(paste(save.path.dir, save.dir, "/padj_hist_", save.dir, ".pdf", sep = ""))
    hist(sleuth_table_wt_filt$qval, breaks = 20, xlab = "adjusted pval by BH", main = save.dir)  
    nbupdown <- length(which(sleuth_table_wt_filt$qval<alpha_thresh))
    nbup <- length(which(sleuth_table_wt_filt$b>0 & sleuth_table_wt_filt$qval<alpha_thresh))
    nbdown <- length(which(sleuth_table_wt_filt$b<0 & sleuth_table_wt_filt$qval<alpha_thresh))
    mtext(paste("nb diff=", nbupdown, "(up=", nbup, ", down=", nbdown, ") on ", dim(sleuth_table_wt_filt)[1], sep = ""), 3)
    abline(v = alpha_thresh, col = "red", lty = 2)     
    legend("topright", legend = paste("pval threshold=", alpha_thresh, sep = ""), col = "red", lty = 2, cex =0.75, bty = "n")
    dev.off()

    ##### Equ/Up/downregulated genes along chromosomes
    tmp.bed <- bed.transcripts
    if (str_detect(save.dir, "rt")){
      if ((save.dir == "cut14_wt_rt") == T | (save.dir == "cut14_cdc15_wt_rt") == T | (save.dir == "cut14_wt_rt_only") == T | (save.dir == "cut14_cdc15_wt_rt_only") == T){
        tmp.bed <- rbind(tmp.bed, bed_cut14_based)
      }else{
        tmp.bed <- rbind(tmp.bed, bed_rrp6D_based)
      }
    }

    for (pref in c("up", "down", "equ")) {
      ind <- which(sleuth_table_wt_filt$qval<alpha_thresh & sleuth_table_wt_filt$b<0)
      if (pref == "up"){
	ind <- which(sleuth_table_wt_filt$qval<alpha_thresh & sleuth_table_wt_filt$b>0)
      }
      if(pref == "equ"){
	ind <- which(sleuth_table_wt_filt$qval>alpha_thresh)
      }
      signif <- sleuth_table_wt_filt[ind, ]
      genes <- signif$target_id
      add <- do.call(rbind, sapply(genes, function(x) tmp.bed[which(tmp.bed$V4 == x)[1], c("V1", "V2", "V3", "V6", "deb")], simplify = F))
      info <- cbind(as.data.frame(signif), add) 	
      #b <- sapply(genes, function(x) sleuth_table_wt_filt[which(sleuth_table_wt_filt$target_id == x), ]$b)
      #info <- cbind(genes, gff.transcripts[sapply(genes, function(x) which(keep.names == x)), c("V1", "V4", "V5", "V7", "deb")], b)
      eval(parse(text = paste("info.", pref, ".genes <- info", sep = "")))
      eval(parse(text = paste("chrom.", pref, " <- unique(info$V1)", sep = "")))
    }
    chr.diff.all <- unique(c(chrom.up, chrom.down))
    chr.diff.all <- chr.diff.all[order(chr.diff.all)]
 
    ##### MAplot
    pdf(paste(save.path.dir, save.dir, "/MAplot_diff_expressed_", save.dir, ".pdf", sep = ""))
    plot_ma(so, comp, test_type = 'wt', which_model = "full", sig_level = alpha_thresh, point_alpha = 0.2)
    dev.off()
 
    if (str_detect(save.dir, "_rt_only") == F) {
      ##### Density along chromosomes representation with sliding windows
      png(paste(save.path.dir, save.dir, "/presence_along_chromosomes_", save.dir, ".png", sep = ""), w = 450*length(chr.diff.all), h = 750)# 12)	
      w <- 10000; o <- 1
      den.d <- NULL; pos.d <- NULL; den.nd <- NULL; pos.nd <- NULL; ud <- c(); equ <- c(); up <- c(); dw <- c(); 
      for (chr in chr.diff.all) {
        print(chr)
        tmp.up <- info.up.genes[which(info.up.genes$V1 == chr), ]; pos.genes.up <- unlist(apply(tmp.up, 1, function(x) x["V2"]:x["V3"]))
       tmp.down <- info.down.genes[which(info.down.genes$V1 == chr), ]; pos.genes.down <- unlist(apply(tmp.down, 1, function(x) x["V2"]:x["V3"]))
        tmp.equ <- info.equ.genes[which(info.equ.genes$V1 == chr), ];  pos.genes.not.diff <- unlist(apply(tmp.equ, 1, function(x) x["V2"]:x["V3"]))

        pos.genes.diff <- c(pos.genes.up, pos.genes.down); 
        tab.diff <- ComputeWindowsOverChrom(vect.pos = pos.genes.diff, size.chr = sizes[1, chr]) 
        TS <- zoo(tab.diff)
        roll.diff <- rollapply(TS, width = w, by = o, FUN = mean, align = "right")#, partial = T)
		    
        tab.not.diff <- ComputeWindowsOverChrom(vect.pos = pos.genes.not.diff, size.chr = sizes[1, chr]) 
        TS <- zoo(tab.not.diff)
        roll.not.diff <- rollapply(TS, width = w, by = o, FUN = mean, align = "right")#, partial = T)
		    
        den.d <- c(den.d, roll.diff); pos.d <- c(pos.d, list(as.numeric(names(roll.diff))))
        den.nd <- c(den.nd, roll.not.diff); pos.nd <- c(pos.nd, list(as.numeric(names(roll.not.diff))))
        ud <- c(ud, dim(tmp.up)[1]+dim(tmp.down)[1]); up <- c(up, dim(tmp.up)[1]); dw <- c(dw, dim(tmp.down)[1])
        equ <- c(equ, dim(tmp.equ)[1])
      }
      offset <- sizes[1, chr.diff.all]
      offset <- c(0, cumsum(offset)[1:length(offset)])
      posd <- unlist(lapply(1:length(pos.d), function(x) pos.d[[x]] + offset[x]))
      posnd <- unlist(lapply(1:length(pos.nd), function(x) pos.nd[[x]] + offset[x]))
      lim <- ceiling(max(c(den.d, den.nd)))
      plot(posd, den.d + lim, type = "p", pch = 16, cex = 0.7, main = paste("Presence of features - Chromosomes ", 
        paste(chr.diff.all, collapse = " "), " - ", save.dir, "\n window=", w, "bp - pthresh=", 
        alpha_thresh, sep = ""), xlab = "Position", ylab = paste("Count in window (size=", 
        w, " bp) by overlap of ", o, " bp", sep = ""), xlim = c(w, max(offset)), ylim = c(0, 
        2 * lim), yaxt = "n", col = "purple")
      axis(2, c(0, lim), c(0, lim), col = "black")
      axis(4, c(lim, 2 * lim), c(0, lim), col = "purple", col.ticks = "purple", col.axis = "purple")
      points(posnd, den.nd, type = "p", col = 1, pch = 16, cex = 0.7)
      abline(h = 0, col = "gray")
      abline(h = lim, col = "gray")
      u <- par("usr")
      arrows(u[1], u[3], u[1], u[4], code = 2, xpd = TRUE)
      text(0.95 * offset[2:length(offset)], rep(1.85 * lim, length(offset) - 1), paste(ud, 
        " up (", up, ")\n/downregulated (", dw, ")", sep = ""), col = "darkorchid4", 
        cex = 0.8)
      text(0.95 * offset[2:length(offset)], rep(lim * 0.85, length(offset) - 1), paste(equ, 
        " equally\nexpressed", sep = ""), cex = 0.8, col = "darkgray")
      abline(v = offset, col = "darkgray", lwd = 1.5, lty = 2)
      dev.off()

      ##### Distance between upregulated genes and the next tRNA promoters
      min.dist.up.tRNA <- sapply(1:dim(info.up.genes)[1], function(x) min(abs(bed.tRNA$deb[which(bed.tRNA$V1 == info.up.genes$V1[x])]-info.up.genes$deb[x])))
      min.dist.down.tRNA <- sapply(1:dim(info.down.genes)[1], function(x) min(abs(bed.tRNA$deb[which(bed.tRNA$V1 == info.down.genes$V1[x])]-info.down.genes$deb[x])))
      min.dist.up.down.tRNA <- c(min.dist.down.tRNA, min.dist.up.tRNA)
      min.dist.equal.tRNA <- sapply(1:dim(info.equ.genes)[1], function(x) min(abs(bed.tRNA$deb[which(bed.tRNA$V1 == info.equ.genes$V1[x])]-info.equ.genes$deb[x])))
      df.dist.tRNA <- data.frame(group = c(rep("equally.regulated", length(min.dist.equal.tRNA)), rep("up.down.regulated", length(min.dist.up.down.tRNA))), distance.tRNA = c(min.dist.equal.tRNA, min.dist.up.down.tRNA))
      nb.sample <- length(which(df.dist.tRNA$group == "up.down.regulated"))
      # subsample according to the nb of diff expressed genes
      df.tmp <- rbind(df.dist.tRNA[sample(which(df.dist.tRNA$group == "equally.regulated"), nb.sample, replace = T), ], df.dist.tRNA[which(df.dist.tRNA$group == "up.down.regulated"), ])
      # fit NegBin GLM
      glm.nb.dist.tRNA <- try(glm.nb(distance.tRNA ~ group, data = df.tmp))
      if(class(glm.nb.dist.tRNA) != "try-error") {
        pval.gr <- summary(glm.nb.dist.tRNA )[["coefficients"]][2, 4]	
      }
      pdf(paste(save.path.dir, save.dir, "/dist_next_tRNA_", save.dir, ".pdf", sep = ""))
      ha <- hist(min.dist.equal.tRNA, breaks = 100, plot = F)
      hdiff <- hist(min.dist.up.down.tRNA, breaks = 100, plot = F)
      ylim <- c(0, max(c(hdiff$counts, ha$counts)));  xlim <- c(0, max(c(hdiff$breaks, ha$breaks)))
  
      par(mar=c(5, 4, 4, 6) + 0.1)
      hist(min.dist.up.down.tRNA, breaks = 100, main = paste("Min distance between diff. expressed genes and tRNA promoters\nbreaks=100 - ", save.dir, sep = ""), cex.main = 0.8, xlab = "Distance in bp", ylim = ylim, xlim = xlim, col = adjustcolor("purple", 0.5))
      hist(min.dist.equal.tRNA, breaks = 100, add = T, col = adjustcolor("black", 0.25))
  
      if(class(glm.nb.dist.tRNA) != "try-error") {
        mu.null <- exp(coefficients(glm.nb.dist.tRNA)[1]); mu.1 <- exp(sum(coefficients(glm.nb.dist.tRNA)));
        xequ <- df.dist.tRNA[which(df.dist.tRNA$group == "equally.regulated"), ]$distance.tRNA; xequ <- xequ[order(xequ)]
        xdiff <- df.dist.tRNA[which(df.dist.tRNA$group == "up.down.regulated"), ]$distance.tRNA; xdiff <- xdiff[order(xdiff)]
        pred.null <- dnbinom(xequ, mu = mu.null, size = glm.nb.dist.tRNA$theta)
        pred.1 <- dnbinom(xdiff, mu = mu.1, size = glm.nb.dist.tRNA$theta)
  
        par(new = T)
        plot(xequ, pred.null, col = "black", lwd = 2, xlab = "", ylab = "", type = "l", xaxt = "n", yaxt = "n")
        points(xdiff, pred.1, col = "purple", lwd = 2, xlab = "", ylab = "", type = "l", xaxt = "n", yaxt = "n")
        axis(side = 4, ylim = c(0, 1))
        mtext(side = 4, line = 3, 'Density')
        if (pval.gr<0.05) {
          txt <- paste("\nAt a risk of 5%, tRNA-gene distances are different in diff\nexpressed genes (pred=", round(mu.1), " bp) compared to\nequally expressed genes (pred=", round(mu.null), " bp)\n(pvalue, Wald test=", 
                   format(pval.gr, scientific = T), ", GLM NegBin)", sep = "")
        }else{
          txt <- paste("\nAt a risk of 5%, we cannot say tRNA-gene distances are different in\ndiff expressed genes (pred=", round(mu.1), " bp) compared to\nequally expressed genes(pred=", round(mu.null), " bp)\n(pvalue, Wald test=", format(pval.gr, scientific = T), 
                   ", GLM NegBin)", sep = "")
        }
        text(1.2*xlim[2]/2, max(c(pred.null, pred.1))/2, txt, cex = 0.7)
      } 
      legend("topright", c("differentially expressed genes", "all but differentially expressed genes"), col = c("purple", "darkgray"), lty = c(-1,-1), pch = c(15,15),bty = "n", cex = 0.6)
      dev.off()
    }
    ##### Heatmap of 100 Best DE genes based on pval adj
    sleuth_significant_wt <- dplyr::filter(sleuth_table_wt_filt, qval < alpha_thresh)
    by_wt <- sleuth_significant_wt %>% group_by(target_id)
    by_wt %>% arrange(rev(desc(qval)))
    hmp <- plot_transcript_heatmap(so, by_wt$target_id[1:min(c(100, dim(by_wt)[1]))], units = "tpm", trans = "log", offset = 1) +ggtitle(paste(save.dir," - log(tpm)", sep = ""))+ theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(vjust = 1, 
    size = 8, hjust = 1))
    pdf(paste(save.path.dir, save.dir,"/Heatmap100DE_", save.dir, ".pdf", sep = ""))
    print(hmp)
    dev.off()

    ##### Distance between samples 
    pdf( paste(save.path.dir, save.dir, "/DistSamples_", save.dir, ".pdf", sep = "")) 
    p <- plot_sample_heatmap(so, use_filtered = TRUE, color_high = "white", color_low = "dodgerblue", x_axis_angle = 50)
    print(p)
    dev.off()

    ##### Dispersion plot
    pdf(paste(save.path.dir, save.dir, "/dispersion_estimates_", save.dir, ".pdf", sep = ""))
    full <- plot_mean_var(so, which_model = "full")+ggtitle(paste(save.dir, " - full", sep = ""))
    null <- plot_mean_var(so, which_model = "reduced")+ggtitle(paste(save.dir," - null", sep = ""))
    multiplot(null, full, cols = 2)
    dev.off()

    ##### Density of counts before and after normalization and filtering
    sleuth_matrix_norm <- sleuth_to_matrix(so, 'obs_norm', 'est_counts')$data
    sleuth_matrix_norm <- sleuth_matrix_norm[which(rownames(sleuth_matrix_norm) %in% sleuth_table_wt_filt$target_id), ]
    count.norm <- melt(sleuth_matrix_norm, variable_name = "Samples"); count.norm$samples <- code[1, as.character(count.norm$Var2)]
    count.norm$replicate <- count.norm$Var2
    sleuth_matrix_raw <- sleuth_to_matrix(so, 'obs_raw', 'est_counts')$data
    sleuth_matrix_raw <- sleuth_matrix_raw[which(rownames(sleuth_matrix_raw) %in% sleuth_table_wt_filt$target_id), ]
    count <- melt(sleuth_matrix_raw, variable_name = "Samples"); count$samples <- code[1, as.character(count$Var2)]
    count$replicate <- count$Var2

    g1 <- ggplot(count, aes(log(value+1), fill = samples, linetype=replicate)) + geom_density(alpha = 0.5) + xlab("") + ylab(expression("density"(log(count + 1))))+ggtitle(paste(save.dir, "- raw counts"))
    g2 <-  ggplot(count.norm, aes(log(value+1), fill = samples, linetype=replicate)) + geom_density(alpha = 0.5) + xlab("") + ylab(expression("density"(log(count + 1))))+ggtitle(paste(save.dir, "- norm counts"))
    pdf(paste(save.path.dir, save.dir, "/counts_", save.dir, ".pdf", sep = ""))
    multiplot(g1, g2, cols = 2)
    dev.off()

    ##### Plot log2(norm counts +1) vs condition per chromosome 
    all.equ <- info.equ.genes
    all.diff <- rbind(info.up.genes, info.down.genes)
    all.equ.diff <- rbind(info.up.genes, info.equ.genes, info.down.genes)

    if (str_detect(save.dir, "_rt_only") == F) {
      for (type in c("equ", "diff", "equ.diff")) {
        print(type)
        eval(parse(text = paste("subset <- all.", type, sep = "")))
        subset.genes <- subset$target_id
        count.norm.comp <- count.norm[which(count.norm$Var1 %in% subset.genes), ]
        count.norm.comp$chromosome <- sapply(count.norm.comp$Var1, function(x) subset[which(subset$target_id == x), ]$V1) # replace all.diff by subset 
        ViolinPlot(input.df.with.count = count.norm.comp, ylabel = "log2(normalized count + 1)", save.pdf = paste(save.path.dir, save.dir, "/log2count_", type, "_expressed_", save.dir, ".pdf", sep = "")) 
  
        cond <- unique(count.norm.comp[, c("samples", "chromosome")])
        mean.count.norm.comp <- do.call(rbind, sapply(as.character(unique(count.norm.comp$Var1)), function(x) do.call(rbind, Filter(Negate(is.null), apply(cond, 1, function(y) {tmp <- which(count.norm.comp$samples == y[1] & count.norm.comp$chromosome == y[2] & count.norm.comp$Var1 == x); if(length(tmp)>0){as.vector(c(x, y, mean(count.norm.comp[tmp, "value"])))}}))), simplify = F))
        colnames(mean.count.norm.comp) <- c("genes", "samples", "chromosome", "value")
        mean.count.norm.comp <- as.data.frame(mean.count.norm.comp)
        mean.count.norm.comp$value <- as.numeric(as.character(mean.count.norm.comp$value))
        mean.count.norm.comp$genes <- as.character(mean.count.norm.comp$genes)
        ViolinPlot(input.df.with.count = mean.count.norm.comp, ylabel = "log2(mean normalized count + 1)", save.pdf = paste(save.path.dir, save.dir, "/mean_log2count_", type, "_expressed_", save.dir, ".pdf", sep = ""), T) 

        p <- ggplot(subset, aes(factor(V1), b, alpha = .3)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + theme(legend.position="none")+ geom_jitter(height = 0, width = 0.2) + ylab("b")+xlab("chromosome")+ggtitle(paste(save.dir,"-", type, "genes"))
        pdf(paste(save.path.dir, save.dir, "/FC_", type, "_expressed_", save.dir,  ".pdf", sep = ""))
        print(p) 
        dev.off()
      }
    }
    ##### Is there an enrichement in DE genes in some chromosomes?
    tab.equ <- table(info.equ.genes$V1)
    count.DE.nonDE.chr <- matrix(c(tab.equ, rep(0, length(tab.equ))), ncol = 3, nrow = 2, byrow = T, dimnames = list(c("EE", "DE"), names(tab.equ)))
    tab.up <- table(info.up.genes$V1) 
    count.DE.nonDE.chr[2, names(tab.up)] <- tab.up
    tab.d <- table(info.down.genes$V1)   
    count.DE.nonDE.chr[2, names(tab.d)] <- count.DE.nonDE.chr[2, names(tab.d)]+tab.d

    plot <- data.frame(chromosome = rep(colnames(count.DE.nonDE.chr), 2), count = c(count.DE.nonDE.chr[1, ], count.DE.nonDE.chr[2, ]), type = c(rep(rownames(count.DE.nonDE.chr), each = dim(count.DE.nonDE.chr)[2])))

    g.fisher <- ggplot(data=plot, aes(x=chromosome, y=count, fill=type)) +
    geom_bar(stat="identity", position=position_dodge())+
    geom_text(aes(label=count), vjust=1.6, color="black", position = position_dodge(0.9), size=3.5)+
    scale_fill_brewer(palette="Paired")+ theme_bw() + ggtitle(paste( save.dir,"\nDE/EE genes per chromosome:\np.fisher = ", formatC(fisher.test(count.DE.nonDE.chr,workspace=2e+07)$p.value, digits = 3), sep = ""))

    pdf(paste(save.path.dir, save.dir, "/DE_EE_per_chromosome_", save.dir, ".pdf", sep = ""))
    print(g.fisher)
    dev.off()
  }
}

##### Summary of significant results on DE analysis and RT detection
res_signif_rt <- NULL
i_rt_only <- which(sapply(names(s2c.all.analysis), function(x) str_detect(x, "_rt_only")) == T)
for (i in i_rt_only){
  save.path.dir <- paste("results/readthrough_analysis/sleuth_analysis/", sep = "")
  save.dir <- names(s2c.all.analysis)[i] 
  res <- read.csv(paste(save.path.dir, save.dir, "/sleuth_shrinked_analysis_diff_expressed_", save.dir, ".txt", sep = ""), sep = "\t", stringsAsFactors = F, h = T)
  res <- res[which(sapply(res$target_id, function(x) str_detect(x, ".rt")) == T),]
  for (pt in c(0.05, 0.01, 0.001)){
    nb.equ.res <-  length(which(res$qval >= pt ))
    for (t in c(0, 0.5, 1, 1.5)) {
      nb.down.res <- length(which(res$qval < pt & log2(exp(res$b)) < -t))
      nb.up.res <- length(which(res$qval < pt & log2(exp(res$b)) > t))
      res_signif_rt <- rbind(res_signif_rt, data.frame(type = "sleuth", test = "RT Wald beta 0", analysis = save.dir, pval_thresh = pt, LFC_thresh = t, nb.down = nb.down.res, nb.up = nb.up.res, nb.equ = nb.equ.res, stringsAsFactors = F))
    }
  }
}
write.table(res_signif_rt, "results/readthrough_analysis/sleuth_analysis/sleuth_signif_readthrough_detection_rt_only.txt", col.names = T, row.names = F, sep = "\t")

res_signif_compare_to_deseq2 <- NULL
i_nort <- which(sapply(names(s2c.all.analysis), function(x) str_detect(x, "_rt")) == F)
for (i in i_nort){
  save.path.dir <- paste("results/readthrough_analysis/sleuth_analysis/", sep = "")
  save.dir <- names(s2c.all.analysis)[i] 
  res <- read.csv(paste(save.path.dir, save.dir, "/sleuth_shrinked_analysis_diff_expressed_", save.dir, ".txt", sep = ""), sep = "\t", stringsAsFactors = F, h = T)
  for (pt in c(0.05, 0.01, 0.001)){
    nb.equ.res <-  length(which(res$qval >= pt ))
    for (t in c(0, 0.5, 1, 1.5)) {
      nb.down.res <- length(which(res$qval < pt & log2(exp(res$b)) < -t))
      nb.up.res <- length(which(res$qval < pt & log2(exp(res$b)) > t))
      res_signif_compare_to_deseq2 <- rbind(res_signif_compare_to_deseq2, data.frame(type = "sleuth", test = "Wald beta 0", analysis = save.dir, pval_thresh = pt, LFC_thresh = t, nb.down = nb.down.res, nb.up = nb.up.res, nb.equ = nb.equ.res, stringsAsFactors = F))
    }
  }
}
write.table(res_signif_compare_to_deseq2, "results/readthrough_analysis/sleuth_analysis/sleuth_signif_analysis_classic.txt", col.names = T, row.names = F, sep = "\t")

##### Intersection of DESeq2 results and Kallisto results at different pvalue and log2FC thresholds
pcomp <- c(0.05, 0.01, 0.001)
tcomp <- c(0.5, 1, 1.5)

res_deseq2_cut14 <- read.csv("./results/DESeq2_analysis/2017_12_31_cut14-208_vs_wt/deseq2_shrinked_analysis_diff_expressed_thresholdLFC_0_2017_12_31_cut14-208_vs_wt.txt", sep = "\t", stringsAsFactors = F, h = T)
res_kallisto_cut14 <- read.csv("results/readthrough_analysis/sleuth_analysis/cut14_wt/sleuth_shrinked_analysis_diff_expressed_cut14_wt.txt", sep = "\t", stringsAsFactors = F, h = T)

res_deseq2_rrp6D <- read.csv("./results/DESeq2_analysis/2017_12_31_rrp6D_vs_wt/deseq2_shrinked_analysis_diff_expressed_thresholdLFC_0_2017_12_31_rrp6D_vs_wt.txt", sep = "\t", stringsAsFactors = F, h = T)
res_kallisto_rrp6D <- read.csv("results/readthrough_analysis/sleuth_analysis/rrp6D_wt/sleuth_shrinked_analysis_diff_expressed_rrp6D_wt.txt", sep = "\t", stringsAsFactors = F, h = T) 

res_deseq2_cut14_cdc15 <- read.csv("./results/DESeq2_analysis/2017_12_31_cut14-208_cdc15-118_vs_wt/deseq2_shrinked_analysis_diff_expressed_thresholdLFC_0_2017_12_31_cut14-208_cdc15-118_vs_wt.txt", sep = "\t", stringsAsFactors = F, h = T) 
res_kallisto_cut14_cdc15 <- read.csv("results/readthrough_analysis/sleuth_analysis/cut14_cdc15_wt/sleuth_shrinked_analysis_diff_expressed_cut14_cdc15_wt.txt", sep = "\t", stringsAsFactors = F, h = T)

res_deseq2_cdc15 <- read.csv("./results/DESeq2_analysis/2017_12_31_cdc15-118_vs_wt/deseq2_shrinked_analysis_diff_expressed_thresholdLFC_0_2017_12_31_cdc15-118_vs_wt.txt", sep = "\t", stringsAsFactors = F, h = T)
res_kallisto_cdc15 <- read.csv("results/readthrough_analysis/sleuth_analysis/cdc15_wt/sleuth_shrinked_analysis_diff_expressed_cdc15_wt.txt", sep = "\t", stringsAsFactors = F, h = T) 

for (suff in c("cut14", "rrp6D", "cut14_cdc15", "cdc15")){
  for (p in pcomp){
    for (t in tcomp) {
      eval(parse(text = paste("kall <- res_kallisto_", suff, "$target_id[which(res_kallisto_", suff, "$qval < ", p," & res_kallisto_", suff, "$log2FC>", t, ")]", sep = "")))
      eval(parse(text = paste("deseq <- rownames(res_deseq2_", suff, ")[which(res_deseq2_", suff, "$padj<", p, " & res_deseq2_", suff, "$log2FoldChange>", t, ")]", sep = "")))
      inter <- length(intersect(sapply(deseq, function(x) strsplit(x, ":")[[1]][2]),
sapply(kall, function(x) paste(head(strsplit(strsplit(x, ":")[[1]][2], ".", fixed = T)[[1]], -1), collapse = "."))))
      pdf(paste("results/readthrough_analysis/sleuth_analysis/", suff, "_pthresh",p, "_log2FCthresh",t, ".pdf", sep = ""))
      draw.pairwise.venn(area1 = length(deseq), area2 = length(kall), cross.area = inter, category = c(paste(suff, "\nDESeq2 up\n(p=", p, ", log2FC=", t, ")", sep = ""), paste(suff, "\nKallisto up\n(p=", p, ", log2FC=", t, ")", sep = "")))
      dev.off()
    }
  }
}


##### Compute average coverage per sample
cov_cut14 <- c(); cov_wt <- c(); cov_cut14_cdc15 <- c(); cov_rrp6D <- c();
for (i in 1:3){
  tmp <- read.csv(paste("results/mapping/mapped/cut14-208/cut14-208_R", i, "_trim_sort_idxstats_report.txt", sep = ""),sep = "\t", h = F)
  cov_cut14 <- c(cov_cut14, sum(tmp$V2*50)/sum(tmp$V3))

  tmp <- read.csv(paste("results/mapping/mapped/wt/wt_R", i, "_trim_sort_idxstats_report.txt", sep = ""),sep = "\t", h = F)
  cov_wt <- c(cov_wt, sum(tmp$V2*50)/sum(tmp$V3))

  tmp <- read.csv(paste("results/mapping/mapped/rrp6D/rrp6D_R", i, "_trim_sort_idxstats_report.txt", sep = ""),sep = "\t", h = F)
  cov_rrp6D <- c(cov_rrp6D, sum(tmp$V2*50)/sum(tmp$V3))

  tmp <- read.csv(paste("results/mapping/mapped/cut14-208_cdc15-118/cdc15_cut14_R", i, "_trim_sort_idxstats_report.txt", sep = ""),sep = "\t", h = F)
  cov_cut14_cdc15 <- c(cov_cut14_cdc15, sum(tmp$V2*50)/sum(tmp$V3))
}


##### Filter putative RT based on log2FC of 1.5 and pvalue of 0.001
pthresh_rt <- 0.001
threshLFC_rt <- 1.5
n <- 30 # nb bins within gene to make the metagene within the gene 
bin <- 100 #
b <- 3000

# read coverage of WT reference tripicates 
ref_R1 <- tbl_df(fread("results/readthrough_analysis/metagene_readthrough/wt_R1.bed", sep = "\t", h = F))
ref_R2 <- tbl_df(fread("results/readthrough_analysis/metagene_readthrough/wt_R2.bed", sep = "\t", h = F))
ref_R3 <- tbl_df(fread("results/readthrough_analysis/metagene_readthrough/wt_R3.bed", sep = "\t", h = F))
chrom <- unique(ref_R1$V1)
for (chr in chrom){
  for (k in 1:3) {
    print(k)
    eval(parse(text = paste("ref_", chr, "_R", k, " <- filter(ref_R", k, ", V1 == \"", chr, "\")", sep = "")))
  }
}
rm(ref_R1); rm(ref_R2); rm(ref_R3)

for (i in i_rt_only){
  name_ref <- names(bed.used)[i]
  ref <- bed_cut14_based
  if(name_ref == "rrp6D_wt"){
    ref <- bed_rrp6D_based
  }
  s2c <- s2c.all.analysis[[i]]

  save.path.dir <- paste("results/readthrough_analysis/sleuth_analysis/", sep = "")
  save.dir <- names(s2c.all.analysis)[i] 
  res <- read.csv(paste(save.path.dir, save.dir, "/sleuth_shrinked_analysis_diff_expressed_", save.dir, ".txt", sep = ""), sep = "\t", stringsAsFactors = F, h = T)
  res.filt <- res[which(res$qval < pthresh_rt & log2(exp(res$b)) > threshLFC_rt), ]
  ref.filt.rt <- res.filt$target_id[which(sapply(res.filt$target_id, function(x) str_detect(x, ".rt")) == T)]
  ref.rt.init <- res$target_id[which(sapply(res$target_id, function(x) str_detect(x, ".rt")) == T)]

  ref.target.rt <- ref[sapply(ref.filt.rt, function(x) which(ref$V4 == x)), ]
  ref.target <- ref[sapply(ref.filt.rt, function(x) which(ref$V4 == x)+1), ]

  # Histogram of readthrough length
  rt.pos <- ref.target.rt[which(ref.target.rt$V6 == "+"), ]
  g.pos <- ref.target[which(ref.target$V6 == "+"), ]
  l.pos <- rt.pos$V3-g.pos$V3+1
  rt.neg <- ref.target.rt[which(ref.target.rt$V6 == "-"), ]
  g.neg <- ref.target[which(ref.target$V6 == "-"), ]
  l.neg <- rt.neg$V2-g.neg$V2+1
  l <- abs(c(l.neg, l.pos))
  pdf(paste(save.path.dir, save.dir, "/hist_readthrough_length_", save.dir, ".pdf", sep = ""))
  hist(l, main = paste(save.dir, "\n", paste(paste(names(summary(l)), summary(l), sep = ": "), collapse = " "), "\non ", length(l), " filt. readthroughs", sep = ""), xlab = "Filtered readthrough length", breaks = 20)
  dev.off()

  mutant <- strsplit(name_ref, "_wt")[[1]][1]
  if(length(which(s2c$condition == "cut14"))>0){mutant <- "cut14"}
  if(length(which(s2c$condition == "cut14_cdc15"))>0){mutant <- "cut14_cdc15"}
  mut_R1 <- tbl_df(fread(paste("results/readthrough_analysis/metagene_readthrough/", mutant, "_R1.bed", sep = ""), sep = "\t", h = F))
  mut_R2 <- tbl_df(fread(paste("results/readthrough_analysis/metagene_readthrough/",  mutant, "_R2.bed", sep = ""), sep = "\t", h = F))
  mut_R3 <- tbl_df(fread(paste("results/readthrough_analysis/metagene_readthrough/", mutant, "_R3.bed", sep = ""), sep = "\t", h = F))
  for (chr in unique(ref.target$V1)){
    for (k in 1:3) {
      eval(parse(text = paste("mut_", chr, "_R", k, " <- filter(mut_R", k, ", V1 == \"", chr, "\")", sep = "")))
    }
  }
  rm(mut_R1); rm(mut_R2); rm(mut_R3)

  cov <- list(); cov_rt <- list()
  for (j in 1:dim(ref.target)[1]) {
    print(j)
    ref.gene <- ref.target[j, ]
    chr <- ref.gene$V1
    ref.gene$V2 <- ref.gene$V2+1
    len_rt <- (ref.target.rt[j, ]$V3-ref.target.rt[j, ]$V2)-(ref.gene$V3-ref.gene$V2)

    bins_g <- round(seq(from = ref.gene$V2, to = ref.gene$V3, length = n+1)); if (tail(bins_g,1)!=ref.gene$V3){bins_g[length(bins_g)] <- tmp.gene$V3};
    length_g <- sapply(1:n, function(x)  bins_g[x+1] - bins_g[x])
    mean_norm_gene_ref <- list(); mean_norm_gene_mut <- list()

    bins_ext <- seq(from = 1, to = b+b/n, by = b/n); n_ext <- length(bins_ext)-1
    length_ext <- sapply(1:n_ext, function(x)  bins_ext[x+1] -bins_ext[x])
    mean_norm_5p_ref <- list(); mean_norm_5p_mut <- list()
    mean_norm_3p_ref <- list(); mean_norm_3p_mut <- list()
    mean_norm_3p_ref_rt <- list(); mean_norm_3p_mut_rt <- list()
    bl <- max(c(b, len_rt))
    if(len_rt<bin){bin_rt <- len_rt}else{bin_rt <- bin}

    for (k in 1:3) {
      for (suff in c("ref", "mut")) {
	if(suff == "ref"){corr <- cov_wt[k]}else{eval(parse(text = paste("corr <- cov_", mutant, "[k]", sep = "")))}
        eval(parse(text = paste("tmp <- filter(", suff, "_", chr, "_R", k, ", V2>=ref.gene$V2-", bl, " & V2<=ref.gene$V3+", bl, ")", sep = "")))
        eval(parse(text = paste("tmp_5p <- filter(tmp, V2>=ref.gene$V2-", b, " & V2<=ref.gene$V2)$V3", sep = "")))
        eval(parse(text = paste("tmp_5p_rt <- filter(tmp, V2>=ref.gene$V2-", len_rt, " & V2<=ref.gene$V2)$V3", sep = "")))
        eval(parse(text = paste("tmp_3p <- filter(tmp, V2>=ref.gene$V3 & V2<=ref.gene$V3+", b, ")$V3", sep = "")))
        eval(parse(text = paste("tmp_3p_rt <- filter(tmp, V2>=ref.gene$V3 & V2<=ref.gene$V3+", len_rt, ")$V3", sep = "")))
        if (ref.gene$V8 == "-") {
          eval(parse(text = paste("tmp$V3 <- rev(tmp$V3)", sep = "")))
          eval(parse(text = paste("tmp_5p <- rev(tmp_3p)", sep = "")))
          eval(parse(text = paste("tmp_3p <- rev(tmp_5p)", sep = "")))
          eval(parse(text = paste("tmp_3p_rt <- rev(tmp_5p_rt)", sep = "")))
        }    
        eval(parse(text = paste("count <- sapply(1:n, function(x) mean(filter(tmp, V2 >= bins_g[x] & V2 < bins_g[x+1])$V3/(corr*length_g[x])))", sep = "")))
        eval(parse(text = paste("mean_norm_gene_", suff, " <- c(mean_norm_gene_", suff, ", list(as.vector(count)))", sep = "")))  
        eval(parse(text = paste("count <- rollapply(zoo(tmp_5p/corr), ", bin, ", mean)", sep = "")))
        #eval(parse(text = paste("count <- sapply(1:n, function(x) mean(tmp_5p[bins_ext[x]:(bins_ext[x+1]-1)]/(corr*length_ext[x])))", sep = "")))
        eval(parse(text = paste("mean_norm_5p_", suff, " <- c(mean_norm_5p_", suff, ", list(as.vector(count)*(1/", bin, ")))", sep = "")))
        eval(parse(text = paste("count <- rollapply(zoo(tmp_3p/corr), ", bin, ", mean)", sep = "")))
        #eval(parse(text = paste("count <- sapply(1:n, function(x) mean(tmp_3p[bins_ext[x]:(bins_ext[x+1]-1)]/(corr*length_ext[x])))", sep = "")))
        eval(parse(text = paste("mean_norm_3p_", suff, " <- c(mean_norm_3p_", suff, ", list(as.vector(count)*(1/", bin, ")))", sep = "")))
        eval(parse(text = paste("count <- rollapply(zoo(tmp_3p_rt/corr), ", bin_rt, ", mean)", sep = "")))
        eval(parse(text = paste("mean_norm_3p_", suff, "_rt <- c(mean_norm_3p_", suff, "_rt, list(as.vector(count)*(1/", bin_rt, ")))", sep = "")))
      }
    }
    for (k in 1:3) {
      for (suff in c("ref", "mut")) {
	eval(parse(text = paste(suff, k, " <- list(c(mean_norm_5p_", suff, "[[k]], mean_norm_gene_", suff, "[[k]], mean_norm_3p_", suff, "[[k]]), c(mean_norm_gene_", suff, "[[k]], mean_norm_3p_", suff, "_rt[[k]]))", sep = "")))
      }
    }
    cov <- c(cov, list(data.frame(x = 1:length(ref1[[1]]), ref1 = ref1[[1]], ref2 = ref2[[1]], ref3 = ref3[[1]], mut1 = mut1[[1]], mut2 = mut2[[1]], mut3 = mut3[[1]])))
    cov_rt <- c(cov_rt, list(data.frame(x = 1:length(ref1[[2]]), ref1 = ref1[[2]], ref2 = ref2[[2]], ref3 = ref3[[2]], mut1 = mut1[[2]], mut2 = mut2[[2]], mut3 = mut3[[2]]))) 
  }
  names(cov) <- ref.target$V4
  names(cov_rt) <- ref.target$V4
  saveRDS(cov, paste("results/readthrough_analysis/metagene_readthrough/", mutant, "_wt_metagene_on_mean_pthresh_rt_", pthresh_rt, "_threshLFC_", threshLFC_rt, ".RData", sep = ""))
  saveRDS(cov_rt, paste("results/readthrough_analysis/metagene_readthrough/", mutant, "_wt_metagene_rt_on_mean_pthresh_rt_", pthresh_rt, "_threshLFC_", threshLFC_rt, ".RData", sep = ""))
}

##### Plot metagene profiles
p_kall <- 0.01; threshlog2FC_kall <- 1
mut <- c("cut14", "cut14_cdc15", "rrp6D")
for (j in seq_along(i_rt_only))){
  i <- i_rt_only[j]
  mutant <- mut[j]
  name_ref <- names(bed.used)[i]
  ref <- bed_cut14_based
  if(name_ref == "rrp6D_wt"){
    ref <- bed_rrp6D_based
  }
  s2c <- s2c.all.analysis[[i]]
  save.path.dir <- paste("results/readthrough_analysis/sleuth_analysis/", sep = "")
  save.dir <- names(s2c.all.analysis)[i] 
  res <- read.csv(paste(save.path.dir, save.dir, "/sleuth_shrinked_analysis_diff_expressed_", save.dir, ".txt", sep = ""), sep = "\t", stringsAsFactors = F, h = T)
  res.filt <- res[which(res$qval < pthresh_rt & log2(exp(res$b)) > threshLFC_rt), ]
  ref.filt.rt <- res.filt$target_id[which(sapply(res.filt$target_id, function(x) str_detect(x, ".rt")) == T)]
  ref.rt.init <- res$target_id[which(sapply(res$target_id, function(x) str_detect(x, ".rt")) == T)]
  ref.target.rt <- ref[sapply(ref.filt.rt, function(x) which(ref$V4 == x)), ]
  ref.target <- ref[sapply(ref.filt.rt, function(x) which(ref$V4 == x)+1), ]

  for (suff in c("", "_rt")) {
    cov_meta <- readRDS(paste("results/readthrough_analysis/metagene_readthrough/", mutant, "_wt_metagene", suff, "_on_mean_pthresh_rt_", pthresh_rt, "_threshLFC_", threshLFC_rt, ".RData", sep = ""))
    for (type in c("", "_without_up")){
      if (type != ""){
        transcr_name <- ref.target$V4
        eval(parse(text = paste("kall <- res_kallisto_", mutant, sep = "")))
        not_up <- setdiff(transcr_name, kall$target_id[which(kall$qval < p_kall & log2(exp(b))>threshlog2FC_kall)])
        cov_tmp <- do.call(rbind, cov_meta[not_up])
      }else{
        cov_tmp <- do.call(rbind, cov_meta)
      }
      # Compute mean over all gene of coverage for each bin per replicate and between replicates
      xp <- unique(cov_tmp$x)
      ref1_tmp <- sapply(xp, function(x) mean(cov_tmp[which(cov_tmp$x == x), ]$ref1))
      ref2_tmp <- sapply(xp, function(x) mean(cov_tmp[which(cov_tmp$x == x), ]$ref2))
      ref3_tmp <- sapply(xp, function(x) mean(cov_tmp[which(cov_tmp$x == x), ]$ref3))
      refp <- sapply(xp, function(x) mean(c(ref1_tmp[x], ref2_tmp[x], ref3_tmp[x])))

      mut1_tmp <- sapply(xp, function(x) mean(cov_tmp[which(cov_tmp$x == x), ]$mut1))
      mut2_tmp <- sapply(xp, function(x) mean(cov_tmp[which(cov_tmp$x == x), ]$mut2))
      mut3_tmp <- sapply(xp, function(x) mean(cov_tmp[which(cov_tmp$x == x), ]$mut3))
      mutp <- sapply(xp, function(x) mean(c(mut1_tmp[x], mut2_tmp[x], mut3_tmp[x])))
      ylim_max <- max(c(refp, ref1_tmp, ref2_tmp, ref3_tmp, mut1_tmp,mut2_tmp,mut3_tmp, mutp))
      ylim_min <- min(c(refp, ref1_tmp, ref2_tmp, ref3_tmp, mut1_tmp,mut2_tmp,mut3_tmp, mutp))

      # Change scale of representation to see the gene
      if(suff != ""){
        xtemp <- (length(ref1_tmp)-n)
        xtemp_med <- seq(from = 1, by = bin, length = n+1)
        x <- c(xtemp_med[1:n], seq(from = tail(xtemp_med, 1), by = 1662/2178, length = length(xp)-length(xtemp_med)+1))
      }else{
        xtemp <- (length(ref1_tmp)-n)/2
        xtemp_med <- seq(from = xtemp, by = bin, length = n+2)
        x <- c(1:xtemp, xtemp_med[2:(length(xtemp_med)-1)], seq(from = tail(xtemp_med, 1), by = 1, length = xtemp))
      }
  
      pdf(paste("results/readthrough_analysis/metagene_readthrough/metagene_readthrough", type, "_", mutant, suff, ".pdf", sep = ""), w = 12, h = 8)
      par(mfrow = c(2,2))
      plot(x, ref1_tmp, type = "l", xaxt = "n", main = paste(mutant, " - TSS/TTS +/-", b, " bp\n(p=",pthresh_rt, " - log2FC=", threshLFC_rt, ")", sep = ""), xlab = "", ylim = c(0, ylim_max), ylab = "normalized count by bin size and mean coverage")
      points(x, ref2_tmp, type = "l", lty = 2)
      points(x, ref3_tmp, type = "l", lty = 3)
      points(x, mut1_tmp, type = "l", col = "red")
      points(x, mut2_tmp, type = "l", lty = 2,  col = "red")
      points(x, mut3_tmp, type = "l", lty = 3,  col = "red")
      if (suff != ""){
        axis(1, c(1,tail(xtemp_med, 1)), c("TSS", "TTS"), cex.axis = 0.75)
        abline(v = tail(xtemp_med, 1))
      }else{ 
        axis(1, c(xtemp, xtemp_med[length(xtemp_med)]), c("TSS", "TTS"), cex.axis = 0.75)
        abline(v = c(xtemp, xtemp_med[length(xtemp_med)]))
      }
 
      plot(x, refp, type = "l", xaxt = "n", main = mutant, ylab = "normalized count by bin size and mean coverage", ylim = c(0, ylim_max), xlab = "")
      points(x, mutp, type = "l", col = "red")
      if (suff != ""){
        axis(1, c(1,tail(xtemp_med, 1)), c("TSS", "TTS"), cex.axis = 0.75)
        abline(v = tail(xtemp_med, 1))
      }else{ 
        axis(1, c(xtemp, xtemp_med[length(xtemp_med)]), c("TSS", "TTS"), cex.axis = 0.75)
        abline(v = c(xtemp, xtemp_med[length(xtemp_med)]))
      }

      plot(x, log2(ref1_tmp+0.5), type = "l", xaxt = "n", main = mutant, ylim = log2(c(ylim_min, ylim_max)+0.5), ylab = "log2(normalized count+0.5)", xlab = "")
      points(x, log2(ref2_tmp+0.5), type = "l", lty = 2)
      points(x, log2(ref3_tmp+0.5), type = "l", lty = 3)
      points(x, log2(mut1_tmp+0.5), type = "l", col = "red")
      points(x, log2(mut2_tmp+0.5), type = "l", lty = 2,  col = "red")
      points(x, log2(mut3_tmp+0.5), type = "l", lty = 3,  col = "red")
      if (suff != ""){
        axis(1, c(1,tail(xtemp_med, 1)), c("TSS", "TTS"), cex.axis = 0.75)
        abline(v = tail(xtemp_med, 1))
      }else{ 
        axis(1, c(xtemp, xtemp_med[length(xtemp_med)]), c("TSS", "TTS"), cex.axis = 0.75)
        abline(v = c(xtemp, xtemp_med[length(xtemp_med)]))
      }

      plot(x, log2(refp+0.5), type = "l", xaxt = "n", main = mutant, ylim = log2(c(ylim_min, ylim_max)+0.5), ylab = "log2(normalized count+0.5)", xlab = "")
      points(x, log2(mutp+0.5), type = "l", col = "red")
      if (suff != ""){
        axis(1, c(1,tail(xtemp_med, 1)), c("TSS", "TTS"), cex.axis = 0.75)
        abline(v = tail(xtemp_med, 1))
      }else{ 
        axis(1, c(xtemp, xtemp_med[length(xtemp_med)]), c("TSS", "TTS"), cex.axis = 0.75)
        abline(v = c(xtemp, xtemp_med[length(xtemp_med)]))
      }
      dev.off()
    }
  }     
}   







