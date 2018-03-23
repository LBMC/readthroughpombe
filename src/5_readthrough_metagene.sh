# script to plot metagene of the readthrough events

bin/nextflow src/readthrough_metagene.nf -c src/pipe/conf/readthrough_docker.config \
  --annotation_forward "results/readthrough/peak_calling/2018_02_28_*_peak_calling_forward.bed" \
  --annotation_reverse "results/readthrough/peak_calling/2018_02_28_*_peak_calling_reverse.bed" \
  --bam "results/readthrough/bams/*/*.bam" \
  --mean_size 363 \
  -resume -w /home/laurent/data/work/ -with-dag results/readthrough/metagene_dag.pdf \
  -with-timeline results/readthrough/metagene_timeline

