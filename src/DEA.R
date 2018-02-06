setwd("~/projects/readthroughpombe")
library("tximport")
library("rhdf5")
library("DESeq2")
library("BiocParallel")
register(BiocParallel::MulticoreParam(4))
library("ggplot2")
results_folder <- "results/readthrough/DEA/"
system(paste0("mkdir -p ", results_folder))
RT_folder <- "results/readthrough/quantification/"

load_files <- function(folder, way) {
  files <- list.files(paste0(folder, way))
  files <- files[grepl(".*tsv", files)]
  files <- paste0(folder, way, "/", files)
  names(files) <- gsub(".*PLBD[0-9]{1,2}_(.*)_sorted.*", "\\1", files, perl = T)
  names(files) <- gsub("-", "_", names(files), perl = T)
  kallisto_results <- tximport::tximport(files, type = "kallisto", txOut = TRUE)
  return(kallisto_results)
}

load_files_FR <- function(folder) {
  kallisto_results_F <- load_files(folder, "forward")
  kallisto_results_R <- load_files(folder, "reverse")
  kallisto_results <- kallisto_results_F
  for (results_type in names(kallisto_results)) {
    kallisto_results[[results_type]] <- rbind(
      kallisto_results_F[[results_type]],
      kallisto_results_R[[results_type]]
    )
    rownames(kallisto_results[[results_type]]) <- gsub(
      "transcript:",
      "",
      rownames(kallisto_results[[results_type]])
    )
  }
  return(kallisto_results)
}

extract_condition <- function(
  kallisto_results, condition_a, condition_b, keep_RT = FALSE) {
  sub_kallisto_results <- kallisto_results
  conditions_select <-
    grepl(
      paste0("^", condition_a, "_R[1-3]$"),
      colnames(kallisto_results$counts)
    ) |
    grepl(
      paste0("^", condition_b, "_R[1-3]$"),
      colnames(kallisto_results$counts)
    )
  for (results_type in names(sub_kallisto_results)[-4]) {
    sub_kallisto_results[[results_type]] <-
      sub_kallisto_results[[results_type]][, conditions_select]
    is_RT <- grepl("RT_", rownames(sub_kallisto_results[[results_type]]))
    if (keep_RT) {
      sub_kallisto_results[[results_type]] <-
        sub_kallisto_results[[results_type]][ is_RT, ]
    } else {
      sub_kallisto_results[[results_type]] <-
        sub_kallisto_results[[results_type]][ !is_RT, ]
    }
  }
  return(sub_kallisto_results)
}

dea_analysis <- function(
  results_folder,
  sub_kallisto_results,
  condition_a,
  condition_b,
  label = "") {
  analysis <- paste0(condition_a, "_vs_", condition_b, label)
  results_folder <- paste0(results_folder, "/", analysis, "/")
  system(paste0("mkdir -p ", results_folder))
  sampletable <- data.frame(
    condition = gsub(
      "(.*)_R[1-3]", "\\1",
      colnames(sub_kallisto_results$counts), perl = T
    )
  )
  rownames(sampletable) <- colnames(sub_kallisto_results)

  dds <- DESeq2::DESeqDataSetFromTximport(
    sub_kallisto_results,
    sampletable,
    ~condition
  )
  dds <- DESeq2::DESeq(dds, parallel = TRUE)
  ntd <- DESeq2::normTransform(dds)
  res <- results(
    dds,
    lfcThreshold = .5,
    altHypothesis = "greaterAbs",
    parallel = TRUE
  )

  pdf(file = paste0(results_folder, analysis, "_DispEsts.pdf"))
  DESeq2::plotDispEsts(dds)
  dev.off()

  pcadata <- DESeq2::plotPCA(ntd, intgroup = c("condition"), returnData = TRUE)
  percentvar <- round(100 * attr(pcadata, "percentVar"))
  ggplot(pcadata, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentvar[1], "% variance")) +
    ylab(paste0("PC2: ", percentvar[2], "% variance")) +
    coord_fixed() +
    theme_bw()
  ggplot2::ggsave(file = paste0(results_folder, analysis, "_PCA.pdf"))

  pdf(file = paste0(results_folder, analysis, "_MA.pdf"))
  DESeq2::plotMA(res, ylim = c(-2.5, 2.5))
  dev.off()

  write.csv(
    as.data.frame(res),
    file = paste0(results_folder, analysis, ".csv")
  )
}

kallisto_results <- load_files_FR(RT_folder)
condition_a <- "wt"
conditions_b <- gsub("(.*)_R[1-3]", "\\1", names(files), perl = T)
conditions_b <- conditions_b[!grepl(condition_a, conditions_b)]
conditions_b <- levels(as.factor(conditions_b))
condition_b <- conditions_b[1]
for (condition_b in conditions_b) {
  for (keep_RT in c(TRUE, FALSE)) {
    sub_kallisto_results <- extract_condition(
      kallisto_results,
      condition_a,
      condition_b,
      keep_RT
    )
    dea_analysis(
      results_folder, sub_kallisto_results, condition_a, condition_b,
      label = ifelse(keep_RT, "_RT", "")
    )
  }
}
