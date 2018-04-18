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
results$signif <- ifelse(results$padj_cut14 < 0.05, "cut14-208", results$signif)
results$signif <- ifelse(results$padj_rrp6D < 0.05 & results$padj_cut14 < 0.05, "cut14-208 & rrp6D", results$signif)
results$signif_pos <- results$padj_rrp6D < 0.05 &
  results$padj_cut14 < 0.05 &
  results$log2FoldChange_rrp6D > 0 &
  results$log2FoldChange_cut14 > 0 

for (group_name in c("cut14-208", "rrp6D", "cut14-208 & rrp6D")) {
  p <- ggplot(data = results[results$signif == group_name & !is.na(results$signif), ],
      aes(
        x = log2FoldChange_cut14,
        y = log2FoldChange_rrp6D,
      )) +
    geom_point(
      size = 1,
      data = results[results$signif == "no", ],
      aes(
        x = log2FoldChange_cut14,
        y = log2FoldChange_rrp6D
      ),
      colour = "gray" 
    ) +
    geom_point(
      size = 1,
      data = results[results$signif != "no" & !is.na(results$signif), ],
      aes(
        x = log2FoldChange_cut14,
        y = log2FoldChange_rrp6D,
        color = signif
      )
    ) +
    xlim(-5, 10) +
    ylim(-5, 10) +
    geom_smooth(method = "lm") +
    theme_bw() +
    labs(
      title = "differentially expressed genes in rrp6D and cut14-208",
      xlab = "cut14-208 log2FoldChange",
      ylab = "rrp6D log2FoldChange",
      color = "differentially expressed"
    )
  print(p)
  ggsave(
    filename = paste0(
      "results/readthrough/DEA/cut14_208_vs_rrp6D_lm_",
      group_name,
      ".pdf"
    ), plot = p,
    width = 29.7, height = 21, units = "cm", scale = 2
  )
}
results$signif <- as.factor(results$signif)
results$signif <- relevel(results$signif, "no")
model <- lm(data = results, log2FoldChange_rrp6D ~ log2FoldChange_cut14*signif)
summary(model)
rmodel <- MASS::rlm(data = results, log2FoldChange_rrp6D ~ log2FoldChange_cut14*signif)
summary(model)
model$coefficients
rmodel$coefficients

p <- ggplot() +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(
    size = 1,
    data = results[results$signif == "no", ],
    aes(
      x = log2FoldChange_cut14,
      y = log2FoldChange_rrp6D
    ),
    colour = "gray" 
  ) +
  geom_point(
    size = 1,
    data = results[results$signif != "no" & !is.na(results$signif), ],
    aes(
      x = log2FoldChange_cut14,
      y = log2FoldChange_rrp6D,
      color = signif
    )
  ) +
  annotate("text",
    x = 5,
    y = -5,
    label = paste0(
      "italic(cor) == ",
      cor(results$log2FoldChange_rrp6D[!is.na( results$signif )],
      results$log2FoldChange_cut14[!is.na( results$signif )]
      )
    ),
    colour = "black",
    size = 10,
    parse = TRUE
  ) +
  annotate("text",
    x = -3,
    y = 5,
    label = paste0(
      "italic(cor) == ",
      cor(results$log2FoldChange_rrp6D[results$signif %in% "no"],
      results$log2FoldChange_cut14[results$signif %in% "no"]
      )
    ),
    colour = "gray",
    size = 10,
    parse = TRUE
  ) +
  annotate("text",
    x = 5,
    y = ( rmodel$coefficients[2] + rmodel$coefficients[7] ) * 5 + rmodel$coefficients[1] + rmodel$coefficients[4] + 1,
    label = paste0(
      "italic(cor) == ",
      cor(results$log2FoldChange_rrp6D[results$signif %in% "cut14-208 & rrp6D"],
      results$log2FoldChange_cut14[results$signif %in% "cut14-208 & rrp6D"]
      )
    ),
    colour = "#00BA38",
    size = 10,
    parse = TRUE
  ) +
  annotate("text",
    x = 5,
    y = ( rmodel$coefficients[2] + rmodel$coefficients[6] ) * 5 + rmodel$coefficients[1] + rmodel$coefficients[3] - 1,
    label = paste0(
      "italic(cor) == ",
      cor(results$log2FoldChange_rrp6D[results$signif %in% "cut14-208"],
      results$log2FoldChange_cut14[results$signif %in% "cut14-208"]
      )
    ),
    colour = "#F8766D",
    size = 10,
    parse = TRUE
  ) +
  annotate("text",
    x = 5,
    y = ( rmodel$coefficients[2] + rmodel$coefficients[8] ) * 5 + rmodel$coefficients[1] + rmodel$coefficients[5] + 2,
    label = paste0(
      "italic(cor) == ",
      cor(results$log2FoldChange_rrp6D[results$signif %in% "rrp6D"],
      results$log2FoldChange_cut14[results$signif %in% "rrp6D"]
      )
    ),
    colour = "#619CFF",
    size = 10,
    parse = TRUE
  ) +
  geom_abline(
    slope = rmodel$coefficients[2] + rmodel$coefficients[6], 
    intercept = rmodel$coefficients[1] + rmodel$coefficients[3],
    col="#F8766D"
  ) +
  geom_abline(
    slope = rmodel$coefficients[2] + rmodel$coefficients[8], 
    intercept = rmodel$coefficients[1] + rmodel$coefficients[5],
    col="#619CFF"
  ) +
  geom_abline(
    slope = rmodel$coefficients[2] + rmodel$coefficients[7], 
    intercept = rmodel$coefficients[1] + rmodel$coefficients[4],
    col="#00BA38"
  ) +
  xlim(-5, 10) +
  ylim(-5, 10) +
  theme_bw() +
  labs(
    title = "differentially expressed genes in rrp6D and cut14-208 (rlm)",
    xlab = "cut14-208 log2FoldChange",
    ylab = "rrp6D log2FoldChange",
    color = "differentially expressed"
  )
print(p)
ggsave(
  filename = paste0(
    "results/readthrough/DEA/cut14_208_vs_rrp6D_rlm_.pdf"
  ), plot = p,
  width = 29.7, height = 21, units = "cm", scale = 2
)

p <- ggplot() +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(
    size = 1,
    data = results[results$signif == "no", ],
    aes(
      x = log2FoldChange_cut14,
      y = log2FoldChange_rrp6D
    ),
    colour = "gray" 
  ) +
  geom_point(
    size = 1,
    data = results[results$signif != "no" & !is.na(results$signif), ],
    aes(
      x = log2FoldChange_cut14,
      y = log2FoldChange_rrp6D,
      color = signif
    )
  ) +
  annotate("text",
    x = 5,
    y = -5,
    label = paste0(
      "italic(cor) == ",
      cor(results$log2FoldChange_rrp6D[!is.na( results$signif )],
      results$log2FoldChange_cut14[!is.na( results$signif )]
      )
    ),
    colour = "black",
    size = 10,
    parse = TRUE
  ) +
  annotate("text",
    x = -3,
    y = 5,
    label = paste0(
      "italic(cor) == ",
      cor(results$log2FoldChange_rrp6D[results$signif %in% "no"],
      results$log2FoldChange_cut14[results$signif %in% "no"]
      )
    ),
    colour = "gray",
    size = 10,
    parse = TRUE
  ) +
  annotate("text",
    x = 5,
    y = ( model$coefficients[2] + model$coefficients[7] ) * 5 + model$coefficients[1] + model$coefficients[4] + 1,
    label = paste0(
      "italic(cor) == ",
      cor(results$log2FoldChange_rrp6D[results$signif %in% "cut14-208 & rrp6D"],
      results$log2FoldChange_cut14[results$signif %in% "cut14-208 & rrp6D"]
      )
    ),
    colour = "#00BA38",
    size = 10,
    parse = TRUE
  ) +
  annotate("text",
    x = 5,
    y = ( model$coefficients[2] + model$coefficients[6] ) * 5 + model$coefficients[1] + model$coefficients[3] - 1,
    label = paste0(
      "italic(cor) == ",
      cor(results$log2FoldChange_rrp6D[results$signif %in% "cut14-208"],
      results$log2FoldChange_cut14[results$signif %in% "cut14-208"]
      )
    ),
    colour = "#F8766D",
    size = 10,
    parse = TRUE
  ) +
  annotate("text",
    x = 5,
    y = ( model$coefficients[2] + model$coefficients[8] ) * 5 + model$coefficients[1] + model$coefficients[5] + 2,
    label = paste0(
      "italic(cor) == ",
      cor(results$log2FoldChange_rrp6D[results$signif %in% "rrp6D"],
      results$log2FoldChange_cut14[results$signif %in% "rrp6D"]
      )
    ),
    colour = "#619CFF",
    size = 10,
    parse = TRUE
  ) +
  geom_abline(
    slope = model$coefficients[2] + model$coefficients[6], 
    intercept = model$coefficients[1] + model$coefficients[3], col="#F8766D") +
  geom_abline(
    slope = model$coefficients[2] + model$coefficients[8], 
    intercept = model$coefficients[1] + model$coefficients[5], col="#619CFF") +
  geom_abline(
    slope = model$coefficients[2] + model$coefficients[7], 
    intercept = model$coefficients[1] + model$coefficients[4], col="#00BA38") +
  xlim(-5, 10) +
  ylim(-5, 10) +
  theme_bw() +
  labs(
    title = "differentially expressed genes in rrp6D and cut14-208 (lm)",
    xlab = "cut14-208 log2FoldChange",
    ylab = "rrp6D log2FoldChange",
    color = "differentially expressed"
  )
print(p)
ggsave(
  filename = paste0(
    "results/readthrough/DEA/cut14_208_vs_rrp6D_lm_.pdf"
  ), plot = p,
  width = 29.7, height = 21, units = "cm", scale = 2
)
