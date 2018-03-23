setwd("~/projects/readthroughpombe")
library("ggplot2")

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

results <- data.frame(
  gene = results_list[["rrp6D"]][["T"]]$X,
  log2FoldChange_rrp6D = results_list[["rrp6D"]][["T"]]$log2FoldChange,
  log2FoldChange_cut14 = results_list[["cut14_208"]][["T"]]$log2FoldChange,
  stat_rrp6D = results_list[["rrp6D"]][["T"]]$stat,
  stat_cut14 = results_list[["cut14_208"]][["T"]]$stat,
  pvalue_rrp6D = results_list[["rrp6D"]][["T"]]$pvalue,
  pvalue_cut14 = results_list[["cut14_208"]][["T"]]$pvalue,
  padj_rrp6D = results_list[["rrp6D"]][["T"]]$padj,
  padj_cut14 = results_list[["cut14_208"]][["T"]]$padj
)

results$signif <- ifelse(results$padj_rrp6D < 0.05, "rrp6D", "no")
results$signif <- ifelse(results$padj_cut14 < 0.05, "cut14", results$signif)
results$signif <- ifelse(results$padj_rrp6D < 0.05 & results$padj_cut14 < 0.05, "both", results$signif)
results$signif_pos <- results$padj_rrp6D < 0.05 &
  results$padj_cut14 < 0.05 &
  results$log2FoldChange_rrp6D > 0 &
  results$log2FoldChange_cut14 > 0 

model_stat <- lm(stat_cut14 ~ stat_rrp6D + signif:stat_rrp6D, data = results)
summary(model_stat)
model_log <- lm(log2FoldChange_cut14 ~ log2FoldChange_rrp6D + signif_pos + signif_pos:log2FoldChange_rrp6D, data = results)
summary(model_log)

ggplot(data = results, aes(x = stat_cut14, y = stat_rrp6D, color = signif, alpha = signif != "no")) +
  geom_point(size = 1) +
  theme_bw()

ggplot(data = results[results$stat_cut14 != 0 & results$stat_rrp6D != 0, ], aes(x = stat_cut14, y = stat_rrp6D, color = signif, alpha = signif != "no")) +
  geom_point(size = 1) +
  theme_bw()

ggplot(data = results, aes(x = log2FoldChange_cut14, y = log2FoldChange_rrp6D, color = signif, alpha = signif != "no")) +
  geom_point(size = 1) +
  theme_bw()

model <- glm(signif_pos ~ stat_rrp6D*stat_cut14, family = binomial(link='logit'), data = results)
summary(model)
model <- glm(signif_pos ~ stat_rrp6D*stat_cut14, family = binomial(link='logit'), data = results[results$stat_cut14 != 0 & results$stat_rrp6D != 0, ])
summary(model)
model <- glm(signif_pos ~ log2FoldChange_rrp6D*log2FoldChange_cut14, family = binomial(link='logit'), data = results)
summary(model)

results <- data.frame(
  gene = results_list[["rrp6D"]][["RT"]]$X,
  log2FoldChange_rrp6D = results_list[["rrp6D"]][["RT"]]$log2FoldChange,
  log2FoldChange_cut14 = results_list[["cut14_208"]][["RT"]]$log2FoldChange,
  stat_rrp6D = results_list[["rrp6D"]][["RT"]]$stat,
  stat_cut14 = results_list[["cut14_208"]][["RT"]]$stat,
  pvalue_rrp6D = results_list[["rrp6D"]][["RT"]]$pvalue,
  pvalue_cut14 = results_list[["cut14_208"]][["RT"]]$pvalue,
  padj_rrp6D = results_list[["rrp6D"]][["RT"]]$padj,
  padj_cut14 = results_list[["cut14_208"]][["RT"]]$padj
)

results$signif <- ifelse(results$padj_rrp6D < 0.05, "rrp6D", "no")
results$signif <- ifelse(results$padj_cut14 < 0.05, "cut14", results$signif)
results$signif <- ifelse(results$padj_rrp6D < 0.05 & results$padj_cut14 < 0.05, "both", results$signif)
results$signif_pos <- results$padj_rrp6D < 0.05 &
  results$padj_cut14 < 0.05 &
  results$log2FoldChange_rrp6D > 0 &
  results$log2FoldChange_cut14 > 0 

model_stat <- lm(stat_cut14 ~ stat_rrp6D + signif:stat_rrp6D, data = results)
summary(model_stat)
model_log <- lm(log2FoldChange_cut14 ~ log2FoldChange_rrp6D + signif_pos + signif_pos:log2FoldChange_rrp6D, data = results)
summary(model_log)

ggplot(data = results, aes(x = stat_cut14, y = stat_rrp6D, color = signif, alpha = signif != "no")) +
  geom_point(size = 1) +
  theme_bw()

ggplot(data = results[results$stat_cut14 != 0 & results$stat_rrp6D != 0, ], aes(x = stat_cut14, y = stat_rrp6D, color = signif, alpha = signif != "no")) +
  geom_point(size = 1) +
  theme_bw()

ggplot(data = results, aes(x = log2FoldChange_cut14, y = log2FoldChange_rrp6D, color = signif, alpha = signif != "no")) +
  geom_point(size = 1) +
  theme_bw()

model <- glm(signif_pos ~ stat_rrp6D*stat_cut14, family = binomial(link='logit'), data = results)
summary(model)
model <- glm(signif_pos ~ stat_rrp6D*stat_cut14, family = binomial(link='logit'), data = results[results$stat_cut14 != 0 & results$stat_rrp6D != 0, ])
summary(model)
model <- glm(signif_pos ~ log2FoldChange_rrp6D*log2FoldChange_cut14, family = binomial(link='logit'), data = results)
summary(model)
