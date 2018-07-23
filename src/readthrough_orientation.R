rm(list = ls())
library("tidyverse")

bed_colnames <- c("chrom", "chromStart", "chromEnd", "name", "score",
                  "strand", "specie", "type", "score_2", "infos")

find_closest <- function(x, x2, i) {
  chrom2 <- x2 %>% select(chrom) %>% slice(i)
  chromStart2 <- x2 %>% select(chromStart) %>% slice(i) %>% as.integer()
  chromEnd2 <- x2 %>% select(chromEnd) %>% slice(i) %>% as.integer()
  strand2 <- x2 %>% select(strand) %>% slice(i) %>% as.integer()
  pos <- x %>%
    filter(chrom %in% chrom2) %>%
    filter(!(chromStart == chromStart2 & chromEnd <= chromEnd2 & strand == strand2)) %>%
    pull(chromMiddle) %>%
    -chromEnd2 %>%
    abs() %>%
    which.min()
  return(pos)
}

# we want to look if the closest transcript to each transcripts are more in the
# same strand or in the opposite strand
file_T <- "results/readthrough/DEA/T/wt_vs_cut14_208_lfc_greaterAbs_than_0.5/2018_02_15_wt_vs_cut14_208.csv.all.bed"
all_bed <- read_tsv(file_T, col_names = FALSE) %>%
  setNames(bed_colnames) %>%
  mutate(DE = FALSE,
         id = paste0(chrom, chromStart, chromEnd, strand),
         chromMiddle = round(abs(chromStart - chromEnd) / 2) +
         ifelse(strand %in% "+",chromStart, chromEnd),
         closest_pos = NA,
         closest_strand = NA)

for (i in all_bed %>% nrow() %>% seq_len()) {
  all_bed$closest_pos[i] <- all_bed %>%
    find_closest(all_bed, i)
  chrom2 <- all_bed %>% select(chrom) %>% slice(i)
  all_bed$closest_strand[i] <- all_bed %>%
    filter(chrom %in% chrom2) %>%
    slice(all_bed$closest_pos[i]) %>%
    pull(strand)
}

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

RT_cut14_bed <- read_tsv(file_RT, col_names = FALSE) %>%
  setNames(bed_colnames) %>%
  mutate(DE = TRUE,
         id = paste0(chrom, chromStart, chromEnd, strand),
         closest_pos = NA,
         closest_strand = NA)

for (i in RT_cut14_bed %>% nrow() %>% seq_len()) {
  RT_cut14_bed$closest_pos[i] <- all_bed %>%
    find_closest(RT_cut14_bed, i)
  chrom2 <- RT_cut14_bed %>% select(chrom) %>% slice(i)
  RT_cut14_bed$closest_strand[i] <- all_bed %>%
    filter(chrom %in% chrom2) %>%
    slice(RT_cut14_bed$closest_pos[i]) %>%
    pull(strand)
}

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

RT_rrp6D_bed <- read_tsv(file_RT, col_names = FALSE) %>%
  setNames(bed_colnames) %>%
  mutate(DE = TRUE,
         id = paste0(chrom, chromStart, chromEnd, strand),
         closest_pos = NA,
         closest_strand = NA)

for (i in RT_rrp6D_bed %>% nrow() %>% seq_len()) {
  RT_rrp6D_bed$closest_pos[i] <- all_bed %>%
    find_closest(RT_rrp6D_bed, i)
  chrom2 <- RT_rrp6D_bed %>% select(chrom) %>% slice(i)
  RT_rrp6D_bed$closest_strand[i] <- all_bed %>%
    filter(chrom %in% chrom2) %>%
    slice(RT_rrp6D_bed$closest_pos[i]) %>%
    pull(strand)
}

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

DE_bed_cut14 <- read_tsv(file_cut14, col_names = FALSE) %>%
  setNames(bed_colnames) %>%
  mutate(DE = TRUE,
         id = paste0(chrom, chromStart, chromEnd, strand),
         closest_pos = NA,
         closest_strand = NA)
DE_bed_rrp6 <- read_tsv(file_rrp6, col_names = FALSE) %>%
  setNames(bed_colnames) %>%
  mutate(DE = TRUE,
         id = paste0(chrom, chromStart, chromEnd, strand),
         closest_pos = NA,
         closest_strand = NA)

DE_cut14only_bed <- DE_bed_cut14 %>% filter(!(id %in% (DE_bed_rrp6 %>% pull(id))))

for (i in DE_cut14only_bed %>% nrow() %>% seq_len()) {
  chrom2 <- DE_cut14only_bed %>% select(chrom) %>% slice(i)
  DE_cut14only_bed$closest_pos[i] <- all_bed %>%
    find_closest(DE_cut14only_bed, i)
  DE_cut14only_bed$closest_strand[i] <- all_bed %>%
    filter(chrom %in% chrom2) %>%
    slice(DE_cut14only_bed$closest_pos[i]) %>%
    pull(strand)
}

table(DE_cut14only_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      DE_cut14only_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
print()
table(DE_cut14only_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      DE_cut14only_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
  as.matrix() %>%
  chisq.test()
