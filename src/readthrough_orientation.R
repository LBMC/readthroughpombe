rm(list = ls())
library("tidyverse")
options(warn = 1)

find_closest <- function(x, x2, i) {
  chrom2 <- x2 %>% select(chrom) %>% slice(i) %>% as.character()
  chromStart2 <- x2 %>% select(chromStart) %>% slice(i) %>% as.integer()
  chromEnd2 <- x2 %>% select(chromEnd) %>% slice(i) %>% as.integer()
  strand2 <- x2 %>% select(strand) %>% slice(i) %>% as.character()
  pos <- x %>%
    filter(chrom %in% chrom2) %>%
    filter(!(chromStart == chromStart2 & chromEnd <= chromEnd2 & strand %in% strand2)) %>%
    pull(chromMiddle) %>%
    -chromEnd2 %>%
    abs() %>%
    which.min()
  return(pos)
}

find_all_closest_pos <- function(x, x2) {
  for (i in x2 %>% nrow() %>% seq_len()) {
    x2$closest_pos[i] <- all_bed %>%
      find_closest(x2, i)
    chrom2 <- x2 %>% select(chrom) %>% slice(i)
    x2$closest_strand[i] <- all_bed %>%
      filter(chrom %in% chrom2) %>%
      slice(x2$closest_pos[i]) %>%
      pull(strand)
  }
  return(x2)
}

load_bed_table <- function(file) {
  bed_colnames <- c("chrom", "chromStart", "chromEnd", "name", "score",
                    "strand", "specie", "type", "score_2", "infos")
  read_tsv(file, col_names = FALSE) %>%
  setNames(bed_colnames) %>%
  mutate(DE = FALSE,
         id = paste0(chrom, chromStart, chromEnd, strand),
         chromMiddle = round(abs(chromStart - chromEnd) / 2) +
         ifelse(strand %in% "+",chromStart, chromEnd),
         closest_pos = NA,
         closest_strand = NA) %>% return()
}

surival_prob <- function(coefficients){
  coefficients <- as.data.frame(coefficients)
  coefficients$OR <- exp(coefficients[, 1])
  coefficients$proba <- coefficients$OR / (1 + coefficients$OR)
  for (i in 2:nrow(coefficients)) {
    coefficients$OR[i] <- exp(coefficients[1, 1]) * exp(coefficients[i, 1])
    coefficients$proba[i] <- coefficients$OR[i] / (1 + coefficients$OR[i])
  }
  return(coefficients)
}

glm_binom <- function(x1, x2) {
  for (i in x1 %>% nrow() %>% seq_len()) {
    test <- x2 %>%
      filter(chrom %in% ( x1 %>% slice(i) %>% select(chrom) %>% as.character() ) &
             chromStart == ( x1 %>% slice(i) %>% select(chromStart) %>% as.integer() ) &
             chromEnd >= ( x1 %>% slice(i) %>% select(chromEnd) %>% as.integer() ) &
             strand %in% ( x1 %>% slice(i) %>% select(strand) %>% as.character() )
            ) %>%
      nrow()
    if (test >= 1) {
      x1 <- x1 %>% slice(-i)
    }
  }
  data <- rbind(tibble(strand = x1 %>% pull(strand) %>% as.factor(),
                     closest = x1 %>% pull(closest_strand) %>% as.factor(),
                     type = "transcript"),
              tibble(strand = x2 %>% pull(strand) %>% as.factor(),
                     closest = x2 %>% pull(closest_strand) %>% as.factor(),
                     type = "RT")) %>%
  mutate(type = factor(as.factor(type), levels = c("transcript", "RT")),
         convergent = strand == closest,
         convergent = factor(as.factor(convergent), levels = c(FALSE, TRUE)))
  g <- glm(convergent ~ type, data = data, family = binomial) %>% summary()
  g <- surival_prob(g$coefficients)
  return(g)
}

# we want to look if the closest transcript to each transcripts are more in the
# same strand or in the opposite strand
file_T <- "results/readthrough/DEA/T/wt_vs_cut14_208_lfc_greaterAbs_than_0.5/2018_02_15_wt_vs_cut14_208.csv.all.bed"
all_bed <- load_bed_table(file_T)
all_bed <- all_bed %>% find_all_closest_pos(all_bed)
all_table <- table(all_bed %>% pull(strand) %>% paste0("strand") %>% as.factor(),
      all_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor())

# we want to look if the closest transcript to the RT are more in the same
# strand or in the opposite strand for cut14

file_RT <- "results/readthrough/DEA/RT/wt_vs_cut14_208_RT_lfc_greater_than_0/2018_02_15_wt_vs_cut14_208_RT.csv.bed"
RT_cut14_bed <- load_bed_table(file_RT)
RT_cut14_bed <- all_bed %>% find_all_closest_pos(RT_cut14_bed)
RT_cut14_table <- table(RT_cut14_bed %>% pull(strand) %>% paste0("strand") %>% as.factor(),
      RT_cut14_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor())

# we want to look if the closest transcript to the RT are more in the same
# strand or in the opposite strand for rrp6D

file_RT <- "results/readthrough/DEA/RT/wt_vs_rrp6D_RT_lfc_greater_than_0/2018_02_15_wt_vs_rrp6D_RT.csv.bed"
RT_rrp6D_bed <- load_bed_table(file_RT)
RT_rrp6D_bed <- all_bed %>% find_all_closest_pos(RT_rrp6D_bed)
RT_rrp6D_table <- table(RT_rrp6D_bed %>% pull(strand) %>% paste0("strand") %>% as.factor(),
      RT_rrp6D_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor())

# we want to look if the closest transcript to the cut14 specific DEG are more
# in the same strand or in the opposite strand

file_cut14 <- "results/readthrough/DEA/T/wt_vs_cut14_208_lfc_greaterAbs_than_0.5/2018_02_15_wt_vs_cut14_208.csv.bed"
file_rrp6 <- "results/readthrough/DEA/T/wt_vs_rrp6D_lfc_greaterAbs_than_0.5/2018_02_15_wt_vs_rrp6D.csv.bed"

DE_bed_cut14 <- load_bed_table(file_cut14)
DE_bed_rrp6 <- load_bed_table(file_rrp6)

DE_cut14only_bed <- DE_bed_cut14 %>% filter(!(id %in% (DE_bed_rrp6 %>% pull(id))))
DE_cut14only_bed <- all_bed %>% find_all_closest_pos(DE_cut14only_bed)
DE_cut14only_table <- table(DE_cut14only_bed %>% pull(strand) %>% paste0("strand") %>% as.factor(),
      DE_cut14only_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor())

all_table
RT_cut14_table
RT_rrp6D_table
DE_cut14only_table

glm_binom(all_bed, RT_cut14_bed)
glm_binom(all_bed, RT_rrp6D_bed)
glm_binom(all_bed, DE_cut14only_bed)
