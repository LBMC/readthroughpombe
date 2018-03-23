# script to plot metagene of the readthrough events

mv results/mapping/mapping/2018_01_19_PLBD8_cut14_208_cdc15_118_R2_trim_rev_sort.bam results/mapping/mapping/2018_01_19_PLBD8_cut14-208_cdc15-118_R2_trim_rev_sort.bam
mv results/mapping/mapping/2018_01_19_PLBD8_cut14_208_cdc15_118_R2_trim_rev_sort.bam.bai results/mapping/mapping/2018_01_19_PLBD8_cut14-208_cdc15-118_R2_trim_rev_sort.bam.bai
mv results/mapping/mapping/2018_01_19_PLBD8_cut14_208_cdc15_118_R2_trim_rev.bam results/mapping/mapping/2018_01_19_PLBD8_cut14-208_cdc15_118_R2_trim_rev.bam
mv results/mapping/mapping/2018_01_19_PLBD8_cut14_208_cdc15_118_R2_trim_rev.bam.bai results/mapping/mapping/2018_01_19_PLBD8_cut14-208_cdc15_118_R2_trim_rev.bam.bai
mv results/mapping/mapping/2018_01_19_PLBD4_cdc15_118_R1_trim_rev_sort.bam results/mapping/mapping/2018_01_19_PLBD4_cdc15-118_R1_trim_rev_sort.bam
mv results/mapping/mapping/2018_01_19_PLBD4_cdc15_118_R1_trim_rev_sort.bam.bai results/mapping/mapping/2018_01_19_PLBD4_cdc15-118_R1_trim_rev_sort.bam.bai

bin/nextflow src/readthrough_metagene.nf -c src/pipe/conf/readthrough_docker.config \
  --annotation_forward "results/readthrough/peak_calling/2018_01_31_*_RT_forward.bed" \
  --annotation_reverse "results/readthrough/peak_calling/2018_01_31_*_RT_reverse.bed" \
  --bam "results/mapping/mapping/2018_01_19_*_rev_sort{.bam,.bam.bai}" \
  --mean_size 363 \
  -resume -w /home/laurent/data/work/ -with-dag results/readthrough/metagene_dag.pdf -with-timeline results/readthrough/metagene_timeline

