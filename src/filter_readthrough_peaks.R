##### Filter peaks obtained from MUSIC
require(data.table)

all_analysis <- c("output_cut14_wt_forward", "output_cut14_wt_reverse", 
  "output_cut14_wt_forward_all", "output_cut14_wt_reverse_all", "output_mutants_wt_forward", 
  "output_mutants_wt_reverse", "output_mutants_wt_forward_all", "output_mutants_wt_reverse_all", 
  "output_rrp6D_wt_forward", "output_rrp6D_wt_reverse", "output_cdc15_wt_forward", 
  "output_cdc15_wt_reverse", "output_cut14_cdc15_wt_forward", "output_cut14_cdc15_wt_reverse")
threshold <- 0.05

##### Concatenate all potential peaks in one bed file
for (analysis in all_analysis) {
  file <- system(paste("ls results/readthrough_analysis/Output_Music/", 
    analysis, "/ERs_[0123456789]*.bed", sep = ""), intern = T)
  system(paste("cp results/readthrough_analysis/Output_Music/", analysis, 
    "/ERs_[0123456789]*.bed results/readthrough_analysis/Output_Music/", 
    analysis, "/ERs_all.bed", sep = ""))
  system(paste("bash src/date.sh results/readthrough_analysis/Output_Music/", 
    analysis, "/ERs_all.bed", sep = ""))
}

##### Clean found peaks by substracted annotation of known features
system("bash src/substract_2_beds.sh")

##### Load annotations for genes + and genes -
gff <- fread("data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr.gff3", 
  h = F, sep = "\t", stringsAsFactors = F)
gff.pos <- as.data.frame(gff[which(gff$V7 == "+" & gff$V3 != "transcript" & 
  gff$V2 == "PomBase" & gff$V1 %in% c("I", "II", "III")), ])
gff.neg <- as.data.frame(gff[which(gff$V7 == "-" & gff$V3 != "transcript" & 
  gff$V2 == "PomBase" & gff$V1 %in% c("I", "II", "III")), ])

##### Prepare annotation file for gene having readthroughand their
##### readthrough
list_analysis <- c("cut14_wt", "cut14_wt_all", "mutants_wt", "mutants_wt_all", 
  "rrp6D_wt", "cdc15_wt", "cut14_cdc15_wt")
all_analysis_list <- list(c("output_cut14_wt_forward", "output_cut14_wt_reverse"), 
  c("output_cut14_wt_forward_all", "output_cut14_wt_reverse_all"), c("output_mutants_wt_forward", 
    "output_mutants_wt_reverse"), c("output_mutants_wt_forward_all", 
    "output_mutants_wt_reverse_all"), c("output_rrp6D_wt_forward", 
    "output_rrp6D_wt_reverse"), c("output_cdc15_wt_forward", "output_cdc15_wt_reverse"), 
  c("output_cut14_cdc15_wt_forward", "output_cut14_cdc15_wt_reverse"))
suffix_list <- list(c("neg", "pos"), c("neg", "pos"), c("neg", "pos"), 
  c("neg", "pos"), c("neg", "pos"), c("neg", "pos"), c("neg", "pos"))

# Keep peaks only matching exactly with 3' of gene + or 5' of gene -
for (l in seq_along(list_analysis)) {
  all_dist.end <- c()
  res <- NULL
  new.annot <- NULL
  all_analysis <- all_analysis_list[[l]]
  suffix <- suffix_list[[l]]
  for (j in seq_along(all_analysis)) {
    analysis <- all_analysis[j]
    suff <- suffix[j]
    #  Open output from substract_2_beds.sh after having made intersection
    # of concatenated raw MUSIC peaks and gene features
    tmp.readthrough <- read.table(paste("results/readthrough_analysis/Output_Music/", 
      analysis, "/ERs_all_intersect_", suff, ".bed", sep = ""), sep = "\t", 
      h = F, stringsAsFactors = F)
    # Select peaks detected in chromosomes I, II, III only
    tmp.readthrough.chr <- tmp.readthrough[which(tmp.readthrough$V1 %in% 
      c("I", "II", "III")), ]
    gff.current <- gff.pos
    gene <- "gene+"
    if (suff == "neg") {
      gff.current <- gff.neg
      gene <- "gene-"
    }
    dist.end <- c()
    # Per chromosome
    for (i in 1:dim(tmp.readthrough.chr)[1]) {
      tmp <- tmp.readthrough.chr[i, ]
      #  subset gff according to chromosome
      tmp.gff <- gff.current[which(gff.current$V1 == tmp$V1), ]
      if (suff == "neg") {
        delta <- tmp.gff$V4 - tmp$V3
      } else {
        delta <- tmp$V2 - tmp.gff$V5
      }
      ind <- which(delta >= 0)
      if (length(ind) > 0) {
        delta <- delta[ind]
        ind.min <- which(delta == min(delta))
        index.gff <- ind[ind.min]
        # Print cases where delta minimal is reached for 2 features
        if (length(ind.min) > 1) {
          l1 <- length(unique(tmp.gff[index.gff, ]$V4))
          l2 <- length(unique(tmp.gff[index.gff, ]$V5))
          if (l1 != 1 | l2 != 1) {
          print(tmp.gff[index.gff, ])
          print(delta[ind.min[1]])
          print("#############")
          }
        }
        dist.end <- c(dist.end, delta[ind.min[1]])
        
        #  Create annotation for readthrough so called 'transcript:gene'.rt'
        if (delta[ind.min[1]] == 0 & length(ind.min) == 1) {
          likely.feature <- tmp.gff[index.gff, ]
          annot.feature.rt <- CreateRTName(likely.feature$V9)
          feature.rt <- likely.feature
          if (suff == "pos") {
          feature.rt$V5 <- tmp$V3
          } else {
          feature.rt$V4 <- tmp$V2
          }
          feature.rt$V9 <- annot.feature.rt
          new.annot <- rbind(new.annot, cbind(feature.rt, data.frame(type = ".rt")))
          new.annot <- rbind(new.annot, cbind(likely.feature, data.frame(type = ".init")))
        }
      }
      rm(delta)
    }
    all_dist.end <- c(all_dist.end, dist.end)
    res <- rbind(res, data.frame(gene_orientation = gene, analysis = analysis, 
      peaks_in_chr = dim(tmp.readthrough.chr)[1], peaks_dist_pos = length(dist.end), 
      peaks_dist_0 = length(which(dist.end == 0)), peaks_dist_5 = length(which(dist.end <= 
        5)), peaks_dist_25 = length(which(dist.end <= 25)), peaks_dist_50 = length(which(dist.end <= 
        50)), stringsAsFactors = F))
  }
  write.table(new.annot, paste("results/readthrough_analysis/Output_Music/annot_readthrough_", 
    list_analysis[[l]], "_forward_reverse.gff3", sep = ""), sep = "\t", 
    col.names = F, row.names = F, quote = F)
  write.table(res, paste("results/readthrough_analysis/Output_Music/summary_annot_readthrough_", 
    list_analysis[[l]], ".txt", sep = ""), sep = "\t", col.names = T, 
    row.names = F, quote = F)
  system(paste("bash src/date.sh results/readthrough_analysis/Output_Music/annot_readthrough_", 
    list_analysis[[l]], "_forward_reverse.gff3", sep = ""))
  system(paste("bash src/date.sh results/readthrough_analysis/Output_Music/summary_annot_readthrough_", 
    list_analysis[[l]], ".txt", sep = ""))
  
  # Plot distance of gene and their putative readthrough
  pdf(paste("results/readthrough_analysis/Output_Music/dist_TTS_peaks_", 
    list_analysis[[l]], "_forward_reverse.pdf", sep = ""), h = 10, 
    w = 7)
  par(mfrow = c(2, 1))
  hist(all_dist.end, breaks = 1500, main = "Distance between TTS and MUSIC peaks:\ncut14 - wt all genes")
  hist(all_dist.end, breaks = 25000, xlim = c(0, 200), main = "Distance between TTS and MUSIC peaks (zoom):\ncut14 - wt all genes")
  dev.off()
  system(paste("bash src/date.sh results/readthrough_analysis/Output_Music/dist_TTS_peaks_", 
    list_analysis[[l]], "_forward_reverse.pdf", sep = ""))
}

#####  Computation of number of readthrough detected before and after
##### substraction with annotation and after filtering with delta
for (l in seq_along(list_analysis)) {
  all_analysis <- all_analysis_list[[l]]
  off <- 1
  plots <- c()
  for (j in seq_along(all_analysis)) {
    analysis <- all_analysis[j]
    all.peaks <- read.csv(paste("results/readthrough_analysis/Output_Music/", 
      analysis, "/ERs_all.bed", sep = ""), sep = "\t", h = F, stringsAsFactors = F)
    file <- system(paste("ls results/readthrough_analysis/Output_Music/", 
      analysis, "/ERs_all_intersect_*", sep = ""), intern = T)
    all.peaks.inter <- read.csv(file, sep = "\t", h = F, stringsAsFactors = F)
    tab <- table(all.peaks$V1)
    plot.tab <- data.frame(chromosome = c("I", "II", "III"), count = as.vector(tab[c("I", 
      "II", "III")]))
    g.tab <- BarplotCountPerChrom(plot.tab, F, paste("After intersection with annotation ", 
      analysis, "\n# found per chromosome (on I, II, III: ", sum(table(all.peaks$V1)[c("I", 
        "II", "III")]), ")", sep = ""))
    eval(parse(text = paste("p", off, " = g.tab", sep = "")))
    plots <- c(plots, off)
    off <- off + 1
    tab <- table(all.peaks.inter$V1)
    plot.tab <- data.frame(chromosome = c("I", "II", "III"), count = as.vector(tab[c("I", 
      "II", "III")]))
    g.tab <- BarplotCountPerChrom(plot.tab, F, paste("After intersection with annotation ", 
      analysis, "\n# found per chromosome (on I, II, III: ", sum(table(all.peaks.inter$V1)[c("I", 
        "II", "III")]), ")", sep = ""))
    eval(parse(text = paste("p", off, " = g.tab", sep = "")))
    plots <- c(plots, off)
    off <- off + 1
  }
  annot.readthrough <- read.csv(paste("results/readthrough_analysis/Output_Music/annot_readthrough_", 
    list_analysis[[l]], "_forward_reverse.gff3", sep = ""), sep = "\t", 
    h = F, stringsAsFactors = F)
  tmp.rt <- annot.readthrough[which(annot.readthrough$V10 == ".rt"), 
    ]
  tab.neg <- table(tmp.rt[which(tmp.rt$V7 == "-"), ]$V1)
  tab.pos <- table(tmp.rt[which(tmp.rt$V7 == "+"), ]$V1)
  
  plot.tab.neg <- data.frame(chromosome = c("I", "II", "III"), count = as.vector(tab.neg[c("I", 
    "II", "III")]))
  g.tab.neg <- BarplotCountPerChrom(plot.tab.neg, F, paste("After filter based on 5' and 3' (genes +)\n# found per chromosome (on I, II, III: ", 
    sum(tab.neg[c("I", "II", "III")]), ")", sep = ""))
  
  plot.tab.pos <- data.frame(chromosome = c("I", "II", "III"), count = as.vector(tab.pos[c("I", 
    "II", "III")]))
  g.tab.pos <- BarplotCountPerChrom(plot.tab.pos, F, paste("After filter based on 5' and 3' (genes +)\n# found per chromosome (on I, II, III: ", 
    sum(tab.pos[c("I", "II", "III")]), ")", sep = ""))
  
  png(paste("results/readthrough_analysis/Output_Music/", list_analysis[l], 
    "_nb_peaks_detected_raw_after_inter.png", sep = ""), h = 750, w = 900)
  eval(parse(text = paste("m = multiplot(", paste(paste("p", plots, sep = ""), 
    collapse = ","), ",g.tab.neg, g.tab.pos, cols=", length(all_analysis) + 
    1, ")", sep = "")))
  dev.off()
  system(paste("bash src/date.sh results/readthrough_analysis/Output_Music/", 
    list_analysis[l], "_nb_peaks_detected_raw_after_inter.png", sep = ""))
}

#####  Put all features together and prepare gff annotation with
##### readthrough added for final analysis
annot.readthrough.cut14.wt.based <- read.csv("results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_forward_reverse.gff3", 
  sep = "\t", h = F)
annot.readthrough.rrp6D.wt.based <- read.csv("results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_forward_reverse.gff3", 
  sep = "\t", h = F)
selected.features <- read.csv("data/ReferenceGenomes/2017_10_03_selected_features_s_pombe.gff3", 
  sep = "\t", h = F)
# Readthroughs specific of one strain
gene.is.not.in.rrp6D <- setdiff(annot.readthrough.cut14.wt.based$V9, annot.readthrough.rrp6D.wt.based$V9)
gene.is.not.in.cut14 <- setdiff(annot.readthrough.rrp6D.wt.based$V9, annot.readthrough.cut14.wt.based$V9)

# Make cut14-wt based + rrp6D annotation
tmp <- annot.readthrough.rrp6D.wt.based[which(annot.readthrough.rrp6D.wt.based$V9 %in% 
  gene.is.not.in.cut14), ]
annot.readthrough.cut14.wt.based.completed <- rbind(annot.readthrough.cut14.wt.based, 
  tmp)
write.table(annot.readthrough.cut14.wt.based.completed[, 1:(dim(annot.readthrough.cut14.wt.based.completed)[2] - 
  1)], "results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse.gff3", 
  sep = "\t", col.names = F, row.names = F, quote = F)
system("bash src/date.sh results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse.gff3")

gene.is.not.complem <- setdiff(selected.features$V9, annot.readthrough.cut14.wt.based.completed$V9)
tmp <- selected.features[which(selected.features$V9 %in% gene.is.not.complem), 
  ]
write.table(rbind(annot.readthrough.cut14.wt.based.completed[, 1:(dim(annot.readthrough.cut14.wt.based.completed)[2] - 
  1)], tmp), "results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.gff3", 
  sep = "\t", col.names = F, row.names = F, quote = F)
system("bash src/date.sh results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.gff3")

# Make rrp6D-wt based + cut14 annotation
tmp <- annot.readthrough.cut14.wt.based[which(annot.readthrough.cut14.wt.based$V9 %in% 
  gene.is.not.in.rrp6D), ]
annot.readthrough.rrp6D.wt.based.completed <- rbind(annot.readthrough.rrp6D.wt.based, 
  tmp)
write.table(annot.readthrough.rrp6D.wt.based.completed[, 1:(dim(annot.readthrough.rrp6D.wt.based.completed)[2] - 
  1)], "results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse.gff3", 
  sep = "\t", col.names = F, row.names = F, quote = F)
system("bash src/date.sh results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse.gff3")

gene.is.not.complem <- setdiff(selected.features$V9, annot.readthrough.rrp6D.wt.based.completed$V9)
tmp <- selected.features[which(selected.features$V9 %in% gene.is.not.complem), 
  ]
write.table(rbind(annot.readthrough.rrp6D.wt.based.completed[, 1:(dim(annot.readthrough.rrp6D.wt.based.completed)[2] - 
  1)], tmp), "results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse_with_genes.gff3", 
  sep = "\t", col.names = F, row.names = F, quote = F)
system("bash src/date.sh results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse_with_genes.gff3")
