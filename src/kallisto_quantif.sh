#!/bin/bash

##### Extract fasta from bed 
bin/bedtools2/bin/bedtools getfasta -fi data/ReferenceGenomes/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa -bed results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.bed -name -fo results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.fasta &&\

bin/bedtools2/bin/bedtools getfasta -fi data/ReferenceGenomes/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa -bed results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_rrp6D_wt_based_cut14_forward_reverse.bed -name -fo results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_rrp6D_wt_based_cut14_forward_reverse.fasta &&\

bin/bedtools2/bin/bedtools getfasta -fi data/ReferenceGenomes/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa -bed data/ReferenceGenomes/2017_12_30_sorted_schizosaccharomyces_pombe.chr.transcript.bed -name -fo results/readthrough_analysis/Output_Music/annot_all_genes.fasta &&\


##### Make kallisto index
kallisto index -k 31 --make-unique -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.index results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.fasta &> results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse_report.txt &&\

kallisto index -k 31 --make-unique -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_rrp6D_wt_based_cut14_forward_reverse.index results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_rrp6D_wt_based_cut14_forward_reverse.fasta &> results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_rrp6D_wt_based_cut14_forward_reverse_report.txt &&\

kallisto index -k 31 --make-unique -i results/readthrough_analysis/Output_Music/annot_all_genes.index results/readthrough_analysis/Output_Music/annot_all_genes.fasta &> results/readthrough_analysis/Output_Music/annot_all_genes_report.txt &&\


##### Make kallisto quantif

# ----- from all genes annotation ----- #
# wt
mkdir results/kallisto_quantif/wt_all_genes_R{1..3} &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/wt_all_genes_R1 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R1_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/wt_all_genes_R3 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R3_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/wt_all_genes_R2 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R2_trim.fastq.gz &&\

# cut14 
mkdir results/kallisto_quantif/cut14_all_genes_R{1..3} &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/cut14_all_genes_R1 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R1_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/cut14_all_genes_R2 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R2_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/cut14_all_genes_R3 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R3_trim.fastq.gz &&\


# cdc15 
mkdir results/kallisto_quantif/cdc15_all_genes_R{1..3} &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/cdc15_all_genes_R1 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/2017_09_18_PLBD4_cdc15_118_R1_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/cdc15_all_genes_R2 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/2017_09_18_PLBD9_cdc15-118_R2_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/cdc15_all_genes_R3 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/2017_09_19_PLBD14_cdc15-118_R3_trim.fastq.gz &&\


# cut14-cdc15
mkdir results/kallisto_quantif/cut14_cdc15_all_genes_R{1..3} &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/cut14_cdc15_all_genes_R1 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*3_cut14-208_cdc15-118_R1_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/cut14_cdc15_all_genes_R2 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*8_cut14_208_cdc15_118_R2_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/cut14_cdc15_all_genes_R3 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*13_cut14-208_cdc15-118_R3_trim.fastq.gz &&\

# rrp6D
mkdir results/kallisto_quantif/rrp6D_all_genes_R{1..3} &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/rrp6D_all_genes_R1 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R1_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/rrp6D_all_genes_R2 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R2_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_all_genes.index -o results/kallisto_quantif/rrp6D_all_genes_R3 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R3_trim.fastq.gz &&\


# ----- from cut14-wt based + rrp6D peak detection ----- #
# wt
#mkdir results/kallisto_quantif/wt_cut14wt_based_+rrp6D_readthrough_with_other_genes_R{1..3} &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.index -o results/kallisto_quantif/wt_cut14wt_based_+rrp6D_readthrough_with_other_genes_R1 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R1_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.index -o results/kallisto_quantif/wt_cut14wt_based_+rrp6D_readthrough_with_other_genes_R3 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R3_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.index -o results/kallisto_quantif/wt_cut14wt_based_+rrp6D_readthrough_with_other_genes_R2 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R2_trim.fastq.gz &&\

# cut14 
#mkdir results/kallisto_quantif/cut14_cut14wt_based_+rrp6D_readthrough_with_other_genes_R{1..3} &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.index -o results/kallisto_quantif/cut14_cut14wt_based_+rrp6D_readthrough_with_other_genes_R1 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R1_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.index -o results/kallisto_quantif/cut14_cut14wt_based_+rrp6D_readthrough_with_other_genes_R2 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R2_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.index -o results/kallisto_quantif/cut14_cut14wt_based_+rrp6D_readthrough_with_other_genes_R3 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R3_trim.fastq.gz &&\

# cut14-cdc15
#mkdir results/kallisto_quantif/cut14_cdc15_cut14wt_based_+rrp6D_readthrough_with_other_genes_R{1..3} &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.index -o results/kallisto_quantif/cut14_cdc15_cut14wt_based_+rrp6D_readthrough_with_other_genes_R1 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*3_cut14-208_cdc15-118_R1_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.index -o results/kallisto_quantif/cut14_cdc15_cut14wt_based_+rrp6D_readthrough_with_other_genes_R2 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*8_cut14_208_cdc15_118_R2_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_cut14_wt_based_rrp6D_forward_reverse.index -o results/kallisto_quantif/cut14_cdc15_cut14wt_based_+rrp6D_readthrough_with_other_genes_R3 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*13_cut14-208_cdc15-118_R3_trim.fastq.gz &&\


# ----- from rrp6D-wt based + cut14 peak detection ----- #
# wt
#mkdir results/kallisto_quantif/wt_rrp6Dwt_based_+cut14_readthrough_with_other_genes_R{1..3} &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_rrp6D_wt_based_cut14_forward_reverse.index -o results/kallisto_quantif/wt_rrp6Dwt_based_+cut14_readthrough_with_other_genes_R1 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R1_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_rrp6D_wt_based_cut14_forward_reverse.index -o results/kallisto_quantif/wt_rrp6Dwt_based_+cut14_readthrough_with_other_genes_R3 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R3_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_rrp6D_wt_based_cut14_forward_reverse.index -o results/kallisto_quantif/wt_rrp6Dwt_based_+cut14_readthrough_with_other_genes_R2 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R2_trim.fastq.gz &&\

# rrp6D
#mkdir results/kallisto_quantif/rrp6D_rrp6Dwt_based_+cut14_readthrough_with_other_genes_R{1..3} &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_rrp6D_wt_based_cut14_forward_reverse.index -o results/kallisto_quantif/rrp6D_rrp6Dwt_based_+cut14_readthrough_with_other_genes_R1 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R1_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_rrp6D_wt_based_cut14_forward_reverse.index -o results/kallisto_quantif/rrp6D_rrp6Dwt_based_+cut14_readthrough_with_other_genes_R2 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R2_trim.fastq.gz &&\
kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_with_other_genes_rrp6D_wt_based_cut14_forward_reverse.index -o results/kallisto_quantif/rrp6D_rrp6Dwt_based_+cut14_readthrough_with_other_genes_R3 --single -l 363.4 -s 85.53354 --bias --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R3_trim.fastq.gz 


