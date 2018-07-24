rm(list = ls())
library("tidyverse")
options(warn = 2)

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

# we want to look if the closest transcript to each transcripts are more in the
# same strand or in the opposite strand
file_T <- "results/readthrough/DEA/T/wt_vs_cut14_208_lfc_greaterAbs_than_0.5/2018_02_15_wt_vs_cut14_208.csv.all.bed"
all_bed <- load_bed_table(file_T)
all_bed <- all_bed %>% find_all_closest_pos(all_bed)

table(all_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      all_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
  print()
table(all_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      all_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
  as.matrix() %>%
  chisq.test()

# we want to look if the closest transcript to the RT are more in the same
# strand or in the opposite strand for cut14

file_RT <- "results/readthrough/DEA/RT/wt_vs_cut14_208_RT_lfc_greater_than_0/2018_02_15_wt_vs_cut14_208_RT.csv.bed"
RT_cut14_bed <- load_bed_table(file_RT)
RT_cut14_bed <- all_bed %>% find_all_closest_pos(RT_cut14_bed)

table(RT_cut14_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      RT_cut14_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
  print()
table(RT_cut14_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      RT_cut14_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
  as.matrix() %>%
  chisq.test()

# we want to look if the closest transcript to the RT are more in the same
# strand or in the opposite strand for rrp6D

file_RT <- "results/readthrough/DEA/RT/wt_vs_rrp6D_RT_lfc_greater_than_0/2018_02_15_wt_vs_rrp6D_RT.csv.bed"
RT_rrp6D_bed <- load_bed_table(file_RT)
RT_rrp6D_bed <- all_bed %>% find_all_closest_pos(RT_rrp6D_bed)

table(RT_rrp6D_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      RT_rrp6D_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
  print()
table(RT_rrp6D_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      RT_rrp6D_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
  as.matrix() %>%
  chisq.test()

# we want to look if the closest transcript to the cut14 specific DEG are more
# in the same strand or in the opposite strand

file_cut14 <- "results/readthrough/DEA/T/wt_vs_cut14_208_lfc_greaterAbs_than_0.5/2018_02_15_wt_vs_cut14_208.csv.bed"
file_rrp6 <- "results/readthrough/DEA/T/wt_vs_rrp6D_lfc_greaterAbs_than_0.5/2018_02_15_wt_vs_rrp6D.csv.bed"

DE_bed_cut14 <- load_bed_table(file_cut14)
DE_bed_rrp6 <- load_bed_table(file_rrp6)

DE_cut14only_bed <- DE_bed_cut14 %>% filter(!(id %in% (DE_bed_rrp6 %>% pull(id))))
DE_cut14only_bed <- all_bed %>% find_all_closest_pos(DE_cut14only_bed)


table(DE_cut14only_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      DE_cut14only_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
print()
table(DE_cut14only_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      DE_cut14only_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
  as.matrix() %>%
  chisq.test()
