# script to plot metagene of the readthrough events

mkdir -p results/readthrough/bams/metagene/
cd results/readthrough/bams/metagene
ln -sf ../cut14/2018_01_19_PLBD11_wt_R3_sorted.bam .
ln -sf ../cut14/2018_01_19_PLBD12_cut14-208_R3_sorted.bam .
ln -sf ../cut14/2018_01_19_PLBD1_wt_R1_sorted.bam .
ln -sf ../cut14/2018_01_19_PLBD2_cut14-208_R1_sorted.bam .
ln -sf ../cut14/2018_01_19_PLBD6_wt_R2_sorted.bam .
ln -sf ../cut14/2018_01_19_PLBD7_cut14-208_R2_sorted.bam .
ln -sf ../rrp6D/2018_01_19_PLBD10_rrp6D_R2_sorted.bam .
ln -sf ../rrp6D/2018_01_19_PLBD15_rrp6D_R3_sorted.bam .
ln -sf ../rrp6D/2018_01_19_PLBD5_rrp6D_R1_sorted.bam .
cd ~/projects/readthroughpombe/

bin/nextflow src/readthrough_metagene.nf -c src/pipe/conf/readthrough_docker.config \
  --annotation_forward "results/readthrough/peak_calling/2018_02_28_*_RT_forward.bed" \
  --annotation_reverse "results/readthrough/peak_calling/2018_02_28_*_RT_reverse.bed" \
  --bam "results/readthrough/bams/metagene/*.bam" \
  -resume -w /home/laurent/data/work/ -with-dag results/readthrough/metagene_dag.pdf \
  -with-timeline results/readthrough/metagene_timeline

