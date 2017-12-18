#!/bin/bash

# Think to make the mkdir to create directories

# Convert gff to bed 
bash src/convertgff_to_bed.sh results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_forward_reverse.gff3 results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_forward_reverse.bed &&\ 

bash src/convertgff_to_bed.sh results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_all_forward_reverse.gff3 results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_all_forward_reverse.bed &&\ 

bash src/convertgff_to_bed.sh results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.gff3 results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.bed &&\ 

bash src/convertgff_to_bed.sh results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.gff3 results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.bed &&\ 

bash src/convertgff_to_bed.sh results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse.gff3 results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse.bed &&\ 

bash src/convertgff_to_bed.sh results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse.gff3 results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse.bed &&\ 

bash src/convertgff_to_bed.sh results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.gff3 results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.bed &&\ 

bash src/convertgff_to_bed.sh results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse_with_genes.gff3 results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse_with_genes.bed &&\ 

bash src/convertgff_to_bed.sh data/ReferenceGenomes/2017_10_03_selected_features_s_pombe.gff3 results/readthrough_analysis/Output_Music/selected_features_s_pombe.bed &&\ 


# Extract fasta from bed 
bin/bedtools2/bin/bedtools getfasta -fi data/ReferenceGenomes/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa -bed results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_forward_reverse.bed -name -fo results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_forward_reverse.fasta &&\

bin/bedtools2/bin/bedtools getfasta -fi data/ReferenceGenomes/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa -bed results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_all_forward_reverse.bed -name -fo results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_all_forward_reverse.fasta &&\

bin/bedtools2/bin/bedtools getfasta -fi data/ReferenceGenomes/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa -bed results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.bed -name -fo results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.fasta &&\

bin/bedtools2/bin/bedtools getfasta -fi data/ReferenceGenomes/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa -bed results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.bed -name -fo results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.fasta &&\


bin/bedtools2/bin/bedtools getfasta -fi data/ReferenceGenomes/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa -bed results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse.bed -name -fo results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse.fasta &&\

bin/bedtools2/bin/bedtools getfasta -fi data/ReferenceGenomes/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa -bed results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse.bed -name -fo results/readthrough_analysis/Output_Music/annot_readthrough_cut14_cdc15_wt_based_forward_reverse.fasta &&\

bin/bedtools2/bin/bedtools getfasta -fi data/ReferenceGenomes/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa -bed results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse.bed -name -fo results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse.fasta &&\


bin/bedtools2/bin/bedtools getfasta -fi data/ReferenceGenomes/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa -bed results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.bed -name -fo results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.fasta &&\

#bin/bedtools2/bin/bedtools getfasta -fi data/ReferenceGenomes/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa -bed results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.bed -name -fo results/readthrough_analysis/Output_Music/annot_readthrough_cut14_cdc15_wt_based_forward_reverse_with_genes.fasta &&\

bin/bedtools2/bin/bedtools getfasta -fi data/ReferenceGenomes/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa -bed results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse_with_genes.bed -name -fo results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse_with_genes.fasta &&\

##### Make kallisto index
/usr/local/bin/kallisto index -k 31 --make-unique -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_forward_reverse.index results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_forward_reverse.fasta &> 2017_11_09_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome_kallisto_report.txt &&\

/usr/local/bin/kallisto index -k 31 --make-unique -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_all_forward_reverse.index results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_all_forward_reverse.fasta &> 2017_11_09_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome_kallisto_cut14_wt_all_report.txt &&\

/usr/local/bin/kallisto index -k 31 --make-unique -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.fasta &> 2017_11_09_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome_kallisto_mutants_wt_report.txt &&\

/usr/local/bin/kallisto index -k 31 --make-unique -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.fasta &> 2017_11_09_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome_kallisto_mutants_wt_all_report.txt &&\

/usr/local/bin/kallisto index -k 31 --make-unique -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.index results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.fasta &> readthrough_cut14_wt_based_forward_reverse_with_genes_report.txt &&\

/usr/local/bin/kallisto index -k 31 --make-unique -i results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse_with_genes.index results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse_with_genes.fasta &> readthrough_rrp6D_wt_based_forward_reverse_with_genes_report.txt &&\


##### Make kallisto quantif

# ----- from cut14/wt peak detection ----- #

# wt 
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_forward_reverse.index -o results/kallisto_quantif/wt_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_forward_reverse.index -o results/kallisto_quantif/wt_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R3_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_forward_reverse.index -o results/kallisto_quantif/wt_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R2_trim.fastq.gz &&\

# cut14 
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_forward_reverse.index -o results/kallisto_quantif/cut14_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_forward_reverse.index -o results/kallisto_quantif/cut14_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R2_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_forward_reverse.index -o results/kallisto_quantif/cut14_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R3_trim.fastq.gz &&\

# ----- from cut14/wt peak detection on all coverages ----- #

# wt 
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_all_forward_reverse.index -o results/kallisto_quantif/wt_all_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_all_forward_reverse.index -o results/kallisto_quantif/wt_all_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R3_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_all_forward_reverse.index -o results/kallisto_quantif/wt_all_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R2_trim.fastq.gz &&\

# cut14
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_all_forward_reverse.index -o results/kallisto_quantif/cut14_all_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_all_forward_reverse.index -o results/kallisto_quantif/cut14_all_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R2_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_all_forward_reverse.index -o results/kallisto_quantif/cut14_all_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R3_trim.fastq.gz &&\


# ----- from mutants/wt peak detection ----- #

# wt
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/wt_mut_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R1_trim.fastq.gz &&\

/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/wt_mut_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R3_trim.fastq.gz &&\

/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/wt_mut_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R2_trim.fastq.gz &&\

# cut14 
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/cut14_mut_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/cut14_mut_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R2_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/cut14_mut_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R3_trim.fastq.gz &&\
    
# rrp6D
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/rrp6D_mut_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/rrp6D_mut_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R2_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/rrp6D_mut_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R3_trim.fastq.gz &&\

# cdc15
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/cdc15_mut_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*4_cdc15_118_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/cdc15_mut_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*9_cdc15-118_R2_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/cdc15_mut_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*14_cdc15-118_R3_trim.fastq.gz &&\

# cut14-cdc15
/usr/local/bin/kallisto quant -i results/read  through_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/cut14_cdc15_mut_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*3_cut14-208_cdc15-118_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/cut14_cdc15_mut_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*8_cut14_208_cdc15_118_R2_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_forward_reverse.index -o results/kallisto_quantif/cut14_cdc15_mut_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*13_cut14-208_cdc15-118_R3_trim.fastq.gz &&\


# ----- from mutants/wt peak detection  on all coverages ----- #

# wt
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/wt_all_mut_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/wt_all_mut_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R3_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/wt_all_mut_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R2_trim.fastq.gz &&\

# cut14 
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/cut14_all_mut_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/cut14_all_mut_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R2_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/cut14_all_mut_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R3_trim.fastq.gz &&\

# rrp6D
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/rrp6D_all_mut_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/rrp6D_all_mut_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R2_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/rrp6D_all_mut_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R3_trim.fastq.gz &&\

# cdc15
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/cdc15_all_mut_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*4_cdc15_118_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/cdc15_all_mut_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*9_cdc15-118_R2_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/cdc15_all_mut_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*14_cdc15-118_R3_trim.fastq.gz &&\

# cut14-cdc15
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/cut14_cdc15_all_mut_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*3_cut14-208_cdc15-118_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/cut14_cdc15_all_mut_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*8_cut14_208_cdc15_118_R2_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_mutants_wt_all_forward_reverse.index -o results/kallisto_quantif/cut14_cdc15_all_mut_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*13_cut14-208_cdc15-118_R3_trim.fastq.gz &&\


# ----- from cut14-wt based + rrp6D peak detection ----- #
# wt
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/wt_cut14wt_based_+rrp6D_with_genes_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/wt_cut14wt_based_+rrp6D_with_genes_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R3_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/wt_cut14wt_based_+rrp6D_with_genes_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R2_trim.fastq.gz &&\

# cut14 
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/cut14_cut14wt_based_+rrp6D_with_genes_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/cut14_cut14wt_based_+rrp6D_with_genes_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R2_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/cut14_cut14wt_based_+rrp6D_with_genes_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*08_R3_trim.fastq.gz &&\

# cut14-cdc15
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/cut14_cdc15_cut14wt_based_+rrp6D_with_genes_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*3_cut14-208_cdc15-118_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/cut14_cdc15_cut14wt_based_+rrp6D_with_genes_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*8_cut14_208_cdc15_118_R2_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_cut14_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/cut14_cdc15_cut14wt_based_+rrp6D_with_genes_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*13_cut14-208_cdc15-118_R3_trim.fastq.gz &&\


# ----- from rrp6D-wt based + cut14 peak detection ----- #
# wt
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/wt_rrp6Dwt_based_+cut14_with_genes_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/wt_rrp6Dwt_based_+cut14_with_genes_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R3_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/wt_rrp6Dwt_based_+cut14_with_genes_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*wt_R2_trim.fastq.gz &&\

# rrp6D
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/rrp6D_rrp6Dwt_based_+cut14_with_genes_R1 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R1_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/rrp6D_rrp6Dwt_based_+cut14_with_genes_R2 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R2_trim.fastq.gz &&\
/usr/local/bin/kallisto quant -i results/readthrough_analysis/Output_Music/annot_readthrough_rrp6D_wt_based_forward_reverse_with_genes.index -o results/kallisto_quantif/rrp6D_rrp6Dwt_based_+cut14_with_genes_R3 --single -l 363.4 -s 85.53354 --bias --single --bootstrap-samples 100 --rf-stranded --threads=4 results/quality_control/trimming/*D_R3_trim.fastq.gz 
