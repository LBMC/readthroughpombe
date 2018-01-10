# generate meta gene

# https://github.com/lh3/bioawk
# cloned the 31st Dec 2017

bin/bioawk/bioawk -c fastx '{ print $name, length($seq) }' <data/ReferenceGenomes/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa > data/ReferenceGenomes/genome.txt &&\

##### compute coverage at bp level

# cut14
bin/bedtools2/bin/bedtools genomecov -d -ibam results/mapping/mapped/cut14-208/2017_10_26_PLBD2_cut14-208_R1_trim_sort.bam -g data/ReferenceGenomes/genome.txt > results/readthrough_analysis/metagene_readthrough/cut14_R1.bed &&\

bin/bedtools2/bin/bedtools genomecov -d -ibam results/mapping/mapped/cut14-208/2017_10_26_PLBD7_cut14-208_R2_trim_sort.bam -g data/ReferenceGenomes/genome.txt > results/readthrough_analysis/metagene_readthrough/cut14_R2.bed &&\

bin/bedtools2/bin/bedtools genomecov -d -ibam results/mapping/mapped/cut14-208/2017_10_26_PLBD12_cut14-208_R3_trim_sort.bam -g data/ReferenceGenomes/genome.txt > results/readthrough_analysis/metagene_readthrough/cut14_R3.bed &&\

# WT
bin/bedtools2/bin/bedtools genomecov -d -ibam results/mapping/mapped/wt/2017_10_26_PLBD1_wt_R1_trim_sort.bam -g data/ReferenceGenomes/genome.txt > results/readthrough_analysis/metagene_readthrough/wt_R1.bed &&\

bin/bedtools2/bin/bedtools genomecov -d -ibam results/mapping/mapped/wt/2017_10_26_PLBD6_wt_R2_trim_sort.bam -g data/ReferenceGenomes/genome.txt > results/readthrough_analysis/metagene_readthrough/wt_R2.bed &&\

bin/bedtools2/bin/bedtools genomecov -d -ibam results/mapping/mapped/wt/2017_10_26_PLBD11_wt_R3_trim_sort.bam -g data/ReferenceGenomes/genome.txt > results/readthrough_analysis/metagene_readthrough/wt_R3.bed &&\

# rrp6D
bin/bedtools2/bin/bedtools genomecov -d -ibam results/mapping/mapped/rrp6D/2017_10_26_PLBD5_rrp6D_R1_trim_sort.bam -g data/ReferenceGenomes/genome.txt > results/readthrough_analysis/metagene_readthrough/rrp6D_R1.bed &&\

bin/bedtools2/bin/bedtools genomecov -d -ibam results/mapping/mapped/rrp6D/2017_10_26_PLBD10_rrp6D_R2_trim_sort.bam -g data/ReferenceGenomes/genome.txt > results/readthrough_analysis/metagene_readthrough/rrp6D_R2.bed &&\

bin/bedtools2/bin/bedtools genomecov -d -ibam results/mapping/mapped/rrp6D/2017_10_26_PLBD15_rrp6D_R3_trim_sort.bam -g data/ReferenceGenomes/genome.txt > results/readthrough_analysis/metagene_readthrough/rrp6D_R3.bed &&\

# cut14 - cdc15
bin/bedtools2/bin/bedtools genomecov -d -ibam results/mapping/mapped/cut14-208_cdc15-118/2017_10_26_PLBD3_cut14-208_cdc15-118_R1_trim_sort.bam -g data/ReferenceGenomes/genome.txt > results/readthrough_analysis/metagene_readthrough/cut14_cdc15_R1.bed &&\

bin/bedtools2/bin/bedtools genomecov -d -ibam results/mapping/mapped/cut14-208_cdc15-118/2017_10_26_PLBD8_cut14_208_cdc15_118_R2_trim_sort.bam -g data/ReferenceGenomes/genome.txt > results/readthrough_analysis/metagene_readthrough/cut14_cdc15_R2.bed &&\

bin/bedtools2/bin/bedtools genomecov -d -ibam results/mapping/mapped/cut14-208_cdc15-118/2017_10_26_PLBD13_cut14-208_cdc15-118_R3_trim_sort.bam -g data/ReferenceGenomes/genome.txt > results/readthrough_analysis/metagene_readthrough/cut14_cdc15_R3.bed &&



##### Compute average coverage per sample
# cut14
samtools idxstats results/mapping/mapped/cut14-208/2017_10_26_PLBD2_cut14-208_R1_trim_sort.bam > results/mapping/mapped/cut14-208/cut14-208_R1_trim_sort_idxstats_report.txt &&\

samtools idxstats results/mapping/mapped/cut14-208/2017_10_26_PLBD7_cut14-208_R2_trim_sort.bam > results/mapping/mapped/cut14-208/cut14-208_R2_trim_sort_idxstats_report.txt &&\

samtools idxstats results/mapping/mapped/cut14-208/2017_10_26_PLBD12_cut14-208_R3_trim_sort.bam > results/mapping/mapped/cut14-208/cut14-208_R3_trim_sort_idxstats_report.txt &&\

# rrp6D
samtools idxstats results/mapping/mapped/rrp6D/2017_10_26_PLBD5_rrp6D_R1_trim_sort.bam > results/mapping/mapped/rrp6D/rrp6D_R1_trim_sort_idxstats_report.txt &&\

samtools idxstats results/mapping/mapped/rrp6D/2017_10_26_PLBD10_rrp6D_R2_trim_sort.bam > results/mapping/mapped/rrp6D/rrp6D_R2_trim_sort_idxstats_report.txt &&\

samtools idxstats results/mapping/mapped/rrp6D/2017_10_26_PLBD15_rrp6D_R3_trim_sort.bam > results/mapping/mapped/rrp6D/rrp6D_R3_trim_sort_idxstats_report.txt &&\

# wt
samtools idxstats results/mapping/mapped/wt/2017_10_26_PLBD1_wt_R1_trim_sort.bam > results/mapping/mapped/wt/wt_R1_trim_sort_idxstats_report.txt &&\

samtools idxstats results/mapping/mapped/wt/2017_10_26_PLBD6_wt_R2_trim_sort.bam > results/mapping/mapped/wt/wt_R2_trim_sort_idxstats_report.txt &&\

samtools idxstats results/mapping/mapped/wt/2017_10_26_PLBD11_wt_R3_trim_sort.bam > results/mapping/mapped/wt/wt_R3_trim_sort_idxstats_report.txt &&\

# cdc15
samtools idxstats results/mapping/mapped/cdc15-118/2017_10_26_PLBD4_cdc15_118_R1_trim_sort.bam > results/mapping/mapped/cdc15-118/cdc15-118_R1_trim_sort_idxstats_report.txt &&\

samtools idxstats results/mapping/mapped/cdc15-118/2017_10_26_PLBD9_cdc15-118_R2_trim_sort.bam > results/mapping/mapped/cdc15-118/cdc15-118_R2_trim_sort_idxstats_report.txt &&\

samtools idxstats results/mapping/mapped/cdc15-118/2017_10_26_PLBD14_cdc15-118_R3_trim_sort.bam > results/mapping/mapped/cdc15-118/cdc15-118_R3_trim_sort_idxstats_report.txt &&\

# cut14-cdc15
samtools idxstats results/mapping/mapped/cut14-208_cdc15-118/2017_10_26_PLBD3_cut14-208_cdc15-118_R1_trim_sort.bam > results/mapping/mapped/cut14-208_cdc15-118/cdc15_cut14_R1_trim_sort_idxstats_report.txt &&\

samtools idxstats results/mapping/mapped/cut14-208_cdc15-118/2017_10_26_PLBD8_cut14_208_cdc15_118_R2_trim_sort.bam > results/mapping/mapped/cut14-208_cdc15-118/cdc15_cut14_R2_trim_sort_idxstats_report.txt &&\

samtools idxstats results/mapping/mapped/cut14-208_cdc15-118/2017_10_26_PLBD13_cut14-208_cdc15-118_R3_trim_sort.bam > results/mapping/mapped/cut14-208_cdc15-118/cdc15_cut14_R3_trim_sort_idxstats_report.txt 
