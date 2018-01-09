##### Filter peaks obtained from MUSIC
require(data.table)
require(ggplot2)
require(stringr)

all_analysis <- c("output_cut14_wt_forward", "output_cut14_wt_reverse", "output_rrp6D_wt_forward", 
  "output_rrp6D_wt_reverse")
threshold <- 0.05

##### Concatenate all potential peaks in one bed file
for (analysis in all_analysis) {
  file <- system(paste("ls results/readthrough_analysis/Output_Music/", 
    analysis, "/ERs_[0123456789]*.bed", sep = ""), intern = T)
  system(paste("cat results/readthrough_analysis/Output_Music/", analysis, 
    "/ERs_*[0123456789].bed > results/readthrough_analysis/Output_Music/", 
    analysis, "/ERs_all.bed", sep = ""))
}

source("src/format_bed.R")

##### Clean found peaks by substracted annotation of known features
system("bash src/substract_2_beds.sh")

##### Load annotations for genes + and genes -
bed <- fread("data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr.bed", 
  h = F, sep = "\t", stringsAsFactors = F)
bed.pos <- as.data.frame(bed[which(bed$V6 == "+" & bed$V8 == "transcript"), ])
bed.neg <- as.data.frame(bed[which(bed$V6 == "-" & bed$V8 == "transcript"), ])

##### Prepare annotation file for gene having readthroughand their
##### readthrough
list_analysis <- c("cut14_wt", "rrp6D_wt")
all_analysis_list <- list(c("output_cut14_wt_forward", "output_cut14_wt_reverse"), 
  c("output_rrp6D_wt_forward", "output_rrp6D_wt_reverse"))
suffix_list <- list(c("neg", "pos"), c("neg", "pos"))

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
    # Open output from substract_2_beds.sh after having made intersection
    # of concatenated raw MUSIC peaks and gene features
    tmp.readthrough <- read.table(paste("results/readthrough_analysis/Output_Music/", 
      analysis, "/ERs_all_intersect_", suff, ".bed", sep = ""), sep = "\t", 
      h = F, stringsAsFactors = F)
    # Select peaks detected in chromosomes I, II, III only
    tmp.readthrough.chr <- tmp.readthrough[which(tmp.readthrough$V1 %in% 
      c("I", "II", "III")), ]
    bed.current <- bed.pos
    gene <- "gene+"
    if (suff == "neg") {
      bed.current <- bed.neg
      gene <- "gene-"
    }
    dist.end <- c()
    for (i in 1:dim(tmp.readthrough.chr)[1]) {
      tmp <- tmp.readthrough.chr[i, ]
      #  subset bed according to chromosome
      tmp.bed <- bed.current[which(bed.current$V1 == tmp$V1), ]
      if (suff == "neg") {
        delta <- tmp.bed$V2 - tmp$V3
      } else {
        delta <- tmp$V2 - tmp.bed$V3
      }
      ind <- which(delta >= 0)
      if (length(ind) > 0) {
        delta <- delta[ind]
        ind.min <- which(delta == min(delta))
        index.bed <- ind[ind.min]
        # Print cases where delta minimal is reached for 2 features
        if (length(ind.min) > 1) {
          l1 <- length(unique(tmp.bed[index.bed, ]$V2))
          l2 <- length(unique(tmp.bed[index.bed, ]$V3))
          if (l1 != 1 | l2 != 1) {
          print(tmp.bed[index.bed, ])
          print(delta[ind.min[1]])
          print("#############")
          }
        }
        dist.end <- c(dist.end, delta[ind.min[1]])
        
        #  Create annotation for readthrough so called 'transcript:gene'.rt'
        if (delta[ind.min[1]] == 0 & length(ind.min) == 1) {
          likely.feature <- tmp.bed[index.bed, ]
          annot.feature.rt <- CreateRTName(likely.feature$V10)
          feature.rt <- likely.feature
          feature.rt$V10 <- annot.feature.rt
          split <- strsplit(annot.feature.rt, ";")[[1]]
          feature.rt$V4 <- na.omit(sapply(split, function(x) strsplit(x, "ID=")[[1]][2]))
          if (suff == "pos") {
            feature.rt$V3 <- tmp$V3
          } else {
            feature.rt$V2 <- tmp$V2
          }
          new.annot <- rbind(new.annot, feature.rt)
          new.annot <- rbind(new.annot, likely.feature)
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
    list_analysis[[l]], "_forward_reverse.bed", sep = ""), sep = "\t", 
    col.names = F, row.names = F, quote = F)
  write.table(res, paste("results/readthrough_analysis/Output_Music/summary_annot_readthrough_", 
    list_analysis[[l]], ".txt", sep = ""), sep = "\t", col.names = T, 
    row.names = F, quote = F)
  #system(paste("bash src/date.sh results/readthrough_analysis/Output_Music/annot_readthrough_", 
  #  list_analysis[[l]], "_forward_reverse.bed", sep = ""))
  #system(paste("bash src/date.sh results/readthrough_analysis/Output_Music/summary_annot_readthrough_", 
  #  list_analysis[[l]], ".txt", sep = ""))
  
  # Plot distance of gene and their putative readthrough
  pdf(paste("results/readthrough_analysis/Output_Music/dist_TTS_peaks_", 
    list_analysis[[l]], "_forward_reverse.pdf", sep = ""), h = 10, 
    w = 7)
  par(mfrow = c(2, 1))
  hist(all_dist.end, breaks = 1500, main = "Distance between TTS and MUSIC peaks:\ncut14 - wt all genes")
  hist(all_dist.end, breaks = 25000, xlim = c(0, 100), main = "Distance between TTS and MUSIC peaks (zoom):\ncut14 - wt all genes")
  dev.off()
  #system(paste("bash src/date.sh results/readthrough_analysis/Output_Music/dist_TTS_peaks_", 
  #  list_analysis[[l]], "_forward_reverse.pdf", sep = ""))
}

#####  Put all features together and prepare gff annotation with
##### transcript annotation from all genes 
bed.transcript <- as.data.frame(bed[which(bed$V8 == "transcript"), ])

##### readthrough added for final analysis
annot.readthrough.cut14.wt.based <- read.csv("results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_forward_reverse.bed", 
  sep = "\t", h = F)
annot.readthrough.rrp6D.wt.based <- read.csv("results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_forward_reverse.bed", 
  sep = "\t", h = F)

# Readthroughs specific of one strain
gene.is.not.in.rrp6D <- setdiff(annot.readthrough.cut14.wt.based$V10, annot.readthrough.rrp6D.wt.based$V10)
gene.is.not.in.cut14 <- setdiff(annot.readthrough.rrp6D.wt.based$V10, annot.readthrough.cut14.wt.based$V10)

##### Make cut14-wt based + rrp6D annotation
tmp <- annot.readthrough.rrp6D.wt.based[which(annot.readthrough.rrp6D.wt.based$V10 %in% 
  gene.is.not.in.cut14), ]
annot.readthrough.cut14.wt.based.completed <- rbind(annot.readthrough.cut14.wt.based, tmp)
write.table(annot.readthrough.cut14.wt.based.completed, "results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_rrp6D_forward_reverse.bed", 
  sep = "\t", col.names = F, row.names = F, quote = F)
#system("bash src/date.sh results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_rrp6D_forward_reverse.bed")

gene.is.not.complem <- setdiff(bed.transcript$V10, annot.readthrough.cut14.wt.based.completed$V10)
tmp <- bed.transcript[which(bed.transcript$V10 %in% gene.is.not.complem), ]
write.table(rbind(annot.readthrough.cut14.wt.based.completed, tmp), "results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.bed", 
  sep = "\t", col.names = F, row.names = F, quote = F)
#system("bash src/date.sh results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.bed")

##### Make rrp6D-wt based + cut14 annotation
tmp <- annot.readthrough.cut14.wt.based[which(annot.readthrough.cut14.wt.based$V10 %in% 
  gene.is.not.in.rrp6D), ]
annot.readthrough.rrp6D.wt.based.completed <- rbind(annot.readthrough.rrp6D.wt.based, tmp)
write.table(annot.readthrough.rrp6D.wt.based.completed, "results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_cut14_forward_reverse.bed", 
  sep = "\t", col.names = F, row.names = F, quote = F)
#system("bash src/date.sh results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_cut14_forward_reverse.bed")

gene.is.not.complem <- setdiff(bed.transcript$V10, annot.readthrough.rrp6D.wt.based.completed$V10)
tmp <- bed.transcript[which(bed.transcript$V10 %in% gene.is.not.complem), ]
write.table(rbind(annot.readthrough.rrp6D.wt.based.completed, tmp), "results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_rrp6D_wt_based_cut14_forward_reverse.bed", 
  sep = "\t", col.names = F, row.names = F, quote = F)
#system("bash src/date.sh results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_rrp6D_wt_based_cut14_forward_reverse_with_genes.bed")
