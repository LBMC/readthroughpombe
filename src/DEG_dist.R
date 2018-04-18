setwd("~/projects/readthroughpombe")
library("ggplot2")

# load DEA results
results_list <- list()
for (condition in c("cut14_208", "rrp6D")) {
  results_list[[condition]] <- list()
  for (analysis in c("RT", "T")) {
    results_list[[condition]][[analysis]] <- read.csv(
      paste0(
        "results/readthrough/DEA/", analysis,
        "/wt_vs_",condition,
        "_lfc_greaterAbs_than_0.5/wt_vs_", condition, ".csv"
      )
    )
  }
}

# load genes positions
bed_list <- list()
for (condition in c("cut14_208", "rrp6D")) {
  bed_list[[condition]] <- list()
  for (analysis in c("RT", "T")) {
    bed_list[[condition]][[analysis]] <- read.csv(
      paste0(
        "results/readthrough/DEA/", analysis,
        "/wt_vs_",condition,
        "_lfc_greaterAbs_than_0.5/wt_vs_", condition, ".csv.all.bed"
      ),
      sep = '\t',
      header = FALSE
    )
  }
}

# list of DEA results
results <- list()
for (analysis in c("RT", "T")) {
  results[[analysis]] <- data.frame(
    gene = results_list[["rrp6D"]][[analysis]]$X,
    log2FoldChange_rrp6D = results_list[["rrp6D"]][[analysis]]$log2FoldChange,
    log2FoldChange_cut14 = results_list[["cut14_208"]][[analysis]]$log2FoldChange,
    stat_rrp6D = results_list[["rrp6D"]][[analysis]]$stat,
    stat_cut14 = results_list[["cut14_208"]][[analysis]]$stat,
    pvalue_rrp6D = results_list[["rrp6D"]][[analysis]]$pvalue,
    pvalue_cut14 = results_list[["cut14_208"]][[analysis]]$pvalue,
    padj_rrp6D = results_list[["rrp6D"]][[analysis]]$padj,
    padj_cut14 = results_list[["cut14_208"]][[analysis]]$padj
  )
}

# table of genes positions
positions <- data.frame(
  chr = bed_list[[ "rrp6D" ]][[ "T" ]][, 1],
  start = bed_list[[ "rrp6D" ]][[ "T" ]][, 2],
  stop = bed_list[[ "rrp6D" ]][[ "T" ]][, 3],
  gene = gsub(
    "transcript:(.*)",
    "\\1",
    bed_list[[ "rrp6D" ]][[ "T" ]][, 4],
    perl = T
  )
)

centromer_pos <- list(
  I = c(3753687, 3789421),
  II = c(1602264, 1644747),
  III = c(1070904, 1137003)
)

centro_distance <- function(chr, start, stop, centro_pos) {
  centro <- centro_pos[[chr]]
  if ((start > centro[1] & start < centro[2]) |
      (stop > centro[1] & stop < centro[2])) {
    return(0)
  }
  # if left of centromer
  if (start < centro[1]) {
    return(min(c(
      centro[1] - start,
      centro[1] - stop
    )))
  } else { # if right of centromer
    return(min(c(
      start - centro[2],
      stop - centro[2]
    )))
  }
}

for (analysis in c("RT", "T")) {
  results[[analysis]]$chr <- NA
  results[[analysis]]$start <- NA
  results[[analysis]]$stop <- NA
  results[[analysis]]$dist <- NA
  for( i in 1:nrow(results[[analysis]])) {
    chr <- as.vector(positions$chr[
      positions$gene %in% results[[analysis]]$gene[i]
    ])
    start <- positions$start[
      positions$gene %in% results[[analysis]]$gene[i]
    ]
    stop <- positions$stop[
      positions$gene %in% results[[analysis]]$gene[i]
    ]
    if (chr %in% names(centromer_pos)) {
      results[[analysis]]$chr[i] <- chr
      results[[analysis]]$start[i] <- start 
      results[[analysis]]$stop[i] <- stop
      results[[analysis]]$distance[i] <- centro_distance(
        chr,
        start,
        stop,
        centromer_pos
      )
    }
  }
}

for (analysis in c("RT", "T")) {
  results[[analysis]]$signif_rrp6D <- FALSE
  results[[analysis]]$signif_rrp6D[results[[analysis]]$padj_rrp6D < 0.05] <- TRUE 
  results[[analysis]]$signif_cut14 <- FALSE
  results[[analysis]]$signif_cut14[results[[analysis]]$padj_cut14 < 0.05] <- TRUE 
  results[[analysis]]$signif_cut14_rrp6D <- results[[analysis]]$signif_cut14 &
    results[[analysis]]$signif_rrp6D
}

ggplot(data = results[["T"]][!is.na( results[["T"]]$padj_rrp6D ), ],
       aes(x = distance, fill = signif_rrp6D)) +
  geom_histogram() +
  theme_bw()

ggplot(data = results[["T"]][!is.na( results[["T"]]$padj_cut14 ), ],
       aes(x = distance, fill = signif_cut14)) +
  geom_histogram() +
  theme_bw()

ggplot(data = results[["T"]][!is.na( results[["T"]]$padj_cut14 ), ],
       aes(x = distance, fill = signif_cut14_rrp6D)) +
  geom_histogram() +
  theme_bw()

model <- MASS::glm.nb(distance ~ signif_cut14,
            data = results[["T"]][!is.na( results[["T"]]$padj_cut14 ), ],
            )
summary(model)
exp(coef(model))

model <- MASS::glm.nb(distance ~ signif_rrp6D,
            data = results[["T"]][!is.na( results[["T"]]$padj_cut14 ), ],
            )
summary(model)
exp(coef(model))

model <- MASS::glm.nb(distance ~ signif_cut14_rrp6D,
            data = results[["T"]][!is.na( results[["T"]]$padj_cut14 ), ],
            )
summary(model)
exp(coef(model))


ggplot(data = results[["RT"]][!is.na( results[["RT"]]$padj_rrp6D ), ],
       aes(x = distance, fill = signif_rrp6D)) +
  geom_histogram() +
  theme_bw()

ggplot(data = results[["RT"]][!is.na( results[["RT"]]$padj_cut14 ), ],
       aes(x = distance, fill = signif_cut14)) +
  geom_histogram() +
  theme_bw()

ggplot(data = results[["RT"]][!is.na( results[["RT"]]$padj_cut14 ), ],
       aes(x = distance, fill = signif_cut14_rrp6D)) +
  geom_histogram() +
  theme_bw()

model <- MASS::glm.nb(distance ~ signif_cut14,
            data = results[["RT"]][!is.na( results[["RT"]]$padj_cut14 ), ],
            )
summary(model)
exp(coef(model))

model <- MASS::glm.nb(distance ~ signif_rrp6D,
            data = results[["RT"]][!is.na( results[["RT"]]$padj_cut14 ), ],
            )
summary(model)
exp(coef(model))

model <- MASS::glm.nb(distance ~ signif_cut14_rrp6D,
            data = results[["RT"]][!is.na( results[["RT"]]$padj_cut14 ), ],
            )
summary(model)
exp(coef(model))





