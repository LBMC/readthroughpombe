##### Filter peaks obtained from MUSIC
require(data.table)

all_analysis <- c("output_cut14_wt_forward", "output_cut14_wt_reverse", 
  "output_cut14_wt_forward_all", "output_cut14_wt_reverse_all", "output_mutants_wt_forward", 
  "output_mutants_wt_reverse", "output_mutants_wt_forward_all", "output_mutants_wt_reverse_all", 
  "output_rrp6D_wt_forward", "output_rrp6D_wt_reverse", "output_cdc15_wt_forward", 
  "output_cdc15_wt_reverse", "output_cut14_cdc15_wt_forward", "output_cut14_cdc15_wt_reverse")
threshold <- 0.05

# Concatenate all potential peaks in one bed file
for (analysis in all_analysis) {
  file <- system(paste("ls results/readthrough_analysis/Output_Music/", 
    analysis, "/ERs_[0123456789]*.bed", sep = ""), intern = T)
  system(paste("cp results/readthrough_analysis/Output_Music/", analysis, 
    "/ERs_[0123456789]*.bed results/readthrough_analysis/Output_Music/", 
    analysis, "/ERs_all.bed", sep = ""))
  system(paste("bash src/date.sh results/readthrough_analysis/Output_Music/", 
    analysis, "/ERs_all.bed", sep = ""))
}

# Clean found peaks by substracted annotation of known features
system("bash src/substract_2_beds.sh")
