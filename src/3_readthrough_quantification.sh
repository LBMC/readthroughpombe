# script to quantifiy the readthrough events

bin/nextflow src/readthrough_quantif.nf -c src/pipe/conf/readthrough_docker.config \
  --annotation_forward "results/readthrough/peak_calling/2018_01_31_*_RT_forward.bed" \
  --annotation_reverse "results/readthrough/peak_calling/2018_01_31_*_RT_reverse.bed" \
  --fasta "data/2017_09_19_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fasta" \
  --bam "results/mapping/mapping/2018_01_19_*_rev_sort.bam" \
  --mean_size 363.4 \
  --sd_size 85.53354 \
  -resume -w /home/laurent/data/work/
