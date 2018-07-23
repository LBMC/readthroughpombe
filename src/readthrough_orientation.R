library("tidyverse")

# we want to look if the closest transcript to the RT are more in the same
# strand or in the opposite strand for cut14

bed_colnames <- c("chrom", "chromStart", "chromEnd", "name", "score",
                  "strand", "specie", "type", "score_2", "infos")
file_RT <- "results/readthrough/DEA/RT/wt_vs_cut14_208_RT_lfc_greater_than_0/2018_02_15_wt_vs_cut14_208_RT.csv.bed"
file_T <- "results/readthrough/DEA/T/wt_vs_cut14_208_lfc_greaterAbs_than_0.5/2018_02_15_wt_vs_cut14_208.csv.DE.csv.all.bed"

RT_cut14_bed <- read_tsv(file_RT, col_names = FALSE) %>%
  setNames(bed_colnames) %>%
  mutate(DE = TRUE,
         closest_pos = NA,
         closest_strand = NA)
all_bed <- read_tsv(file_T, col_names = FALSE) %>%
  setNames(bed_colnames) %>%
  mutate(DE = FALSE,
         chromMiddle = round(abs(chromStart - chromEnd) / 2) +
         ifelse(strand %in% "+",chromStart, chromEnd))

for (i in RT_cut14_bed %>% nrow() %>% seq_len()) {
  chrom <- RT_cut14_bed %>% select(chrom) %>% slice(i)
  pos <- RT_cut14_bed %>% select(chromEnd) %>% slice(i) %>% as.integer()
  RT_cut14_bed$closest_pos[i] <- all_bed %>%
          filter(chrom == chrom) %>%
          pull(chromMiddle) %>%
          -pos %>%
          abs() %>%
          which.min()
  RT_cut14_bed$closest_strand[i] <- all_bed %>%
    filter(chrom == chrom) %>%
    slice(RT_cut14_bed$closest_pos[i]) %>%
    pull(strand)
}

RT_cut14_bed %>% pull(strand) %>% as.factor() %>% table()
RT_cut14_bed %>% pull(closest_strand) %>% as.factor() %>% table()
table(RT_cut14_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      RT_cut14_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
  print()
table(RT_cut14_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      RT_cut14_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
  as.matrix() %>%
  chisq.test()

# we want to look if the closest transcript to the RT are more in the same
# strand or in the opposite strand for rrp6D

bed_colnames <- c("chrom", "chromStart", "chromEnd", "name", "score",
                  "strand", "specie", "type", "score_2", "infos")
file_RT <- "results/readthrough/DEA/RT/wt_vs_rrp6D_RT_lfc_greater_than_0/2018_02_15_wt_vs_rrp6D_RT.csv.bed"
file_T <- "results/readthrough/DEA/T/wt_vs_rrp6D_lfc_greaterAbs_than_0.5/2018_02_15_wt_vs_rrp6D.csv.DE.csv.all.bed"

RT_rrp6D_bed <- read_tsv(file_RT, col_names = FALSE) %>%
  setNames(bed_colnames) %>%
  mutate(DE = TRUE,
         closest_pos = NA,
         closest_strand = NA)
all_bed <- read_tsv(file_T, col_names = FALSE) %>%
  setNames(bed_colnames) %>%
  mutate(DE = FALSE,
         chromMiddle = round(abs(chromStart - chromEnd) / 2) +
         ifelse(strand %in% "+",chromStart, chromEnd))

for (i in RT_rrp6D_bed %>% nrow() %>% seq_len()) {
  chrom <- RT_rrp6D_bed %>% select(chrom) %>% slice(i)
  pos <- RT_rrp6D_bed %>% select(chromEnd) %>% slice(i) %>% as.integer()
  RT_rrp6D_bed$closest_pos[i] <- all_bed %>%
          filter(chrom == chrom) %>%
          pull(chromMiddle) %>%
          -pos %>%
          abs() %>%
          which.min()
  RT_rrp6D_bed$closest_strand[i] <- all_bed %>%
    filter(chrom == chrom) %>%
    slice(RT_rrp6D_bed$closest_pos[i]) %>%
    pull(strand)
}

DE_bed %>% pull(strand) %>% as.factor() %>% table()
DE_bed %>% pull(closest_strand) %>% as.factor() %>% table()
table(DE_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      DE_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
  print()
table(DE_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      DE_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
  as.matrix() %>%
  chisq.test()


# we want to look if the closest transcript to the cut14 specific DEG are more
# in the same strand or in the opposite strand

file_cut14 <- "results/readthrough/DEA/T/wt_vs_cut14_208_lfc_greaterAbs_than_0.5/2018_02_15_wt_vs_cut14_208.csv.bed"
file_rrp6 <- "results/readthrough/DEA/T/wt_vs_rrp6D_lfc_greaterAbs_than_0.5/2018_02_15_wt_vs_rrp6D.csv.bed"

DE_bed_cut14 <- read_tsv(file_cut14, col_names = FALSE) %>%
  setNames(bed_colnames) %>%
  mutate(DE = TRUE,
         id = paste0(chrom, chromStart, chromEnd),
         closest_pos = NA,
         closest_strand = NA)
DE_bed_rrp6 <- read_tsv(file_rrp6, col_names = FALSE) %>%
  setNames(bed_colnames) %>%
  mutate(DE = TRUE,
         id = paste0(chrom, chromStart, chromEnd),
         closest_pos = NA,
         closest_strand = NA)

DE_cut14only_bed <- DE_cut14only_bed_cut14 %>% filter(!(id %in% (DE_cut14only_bed_rrp6 %>% pull(id))))

for (i in DE_cut14only_bed %>% nrow() %>% seq_len()) {
  chrom <- DE_cut14only_bed %>% select(chrom) %>% slice(i)
  pos <- DE_cut14only_bed %>% select(chromEnd) %>% slice(i) %>% as.integer()
  DE_cut14only_bed$closest_pos[i] <- all_bed %>%
          filter(chrom == chrom) %>%
          pull(chromMiddle) %>%
          -pos %>%
          abs() %>%
          which.min()
  DE_cut14only_bed$closest_strand[i] <- all_bed %>%
    filter(chrom == chrom) %>%
    slice(DE_cut14only_bed$closest_pos[i]) %>%
    pull(strand)
}

DE_cut14only_bed %>% pull(strand) %>% as.factor() %>% table()
DE_cut14only_bed %>% pull(closest_strand) %>% as.factor() %>% table()
table(DE_cut14only_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      DE_cut14only_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
print()

table(DE_cut14only_bed %>% pull(strand) %>% paste0("RT") %>% as.factor(),
      DE_cut14only_bed %>% pull(closest_strand) %>% paste0("closest") %>% as.factor()) %>%
  as.matrix() %>%
  chisq.test()
