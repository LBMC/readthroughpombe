# script to plot metagene of the readthrough events

mkdir -p results/readthrough/bams/metagene/
rm results/readthrough/bams/metagene/*.bam*
cd results/readthrough/bams/metagene
find ../../../mapping/bams/ -name "*.bam" | wc -l
find ../../../mapping/bams/ -name "*.bam" | perl -pe "s/(.*(PLBD\d+)-(.+)-Galaxy.*)/\2 \3 \1/g" | sort -k1 > target.txt
find ../../../mapping/mapping/ -name "*.bam" | perl -pe "s/.*(PLBD\d+)_(.*)_R\d.*/\1 \1_\2/g" | sort -k1 > link.txt
join target.txt link.txt | awk '{system("ln -s "$3" "$4"_"$2".bam")}'
ml file_handle
file_handle.py -f *.bam
rm *cdc15*
cd ~/projects/readthroughpombe/

mkdir results/readthrough/metagene/bed/
cd results/readthrough/metagene/bed/
ln -s ../../DEA/RT/wt_vs_cut14_208_RT_lfc_greater_than_0/2018_02_15_wt_vs_cut14_208_RT.csv.bed RT_cut14_208.bed
ln -s ../../DEA/RT/wt_vs_rrp6D_RT_lfc_greater_than_0/2018_02_15_wt_vs_rrp6D_RT.csv.bed RT_rrp6D.bed
file_handle.py -f *.bed
cd ~/projects/readthroughpombe/

bin/nextflow src/readthrough_metagene.nf -c src/pipe/conf/readthrough_docker.config \
  --annotation "results/readthrough/peak_calling/2018_02_28_*_RT_forward.bed" \
  --bam "results/readthrough/bams/metagene/*.bam" \
  -resume -w /home/laurent/data/work/ -with-dag results/readthrough/metagene_dag.pdf \
  -with-timeline results/readthrough/metagene_timeline

ml nextflow
cd /scratch/lmodolo/readthroughpombe/
nextflow src/readthrough_metagene.nf -c src/pipe/conf/readthrough_sge.config \
  --annotation "/scratch/lmodolo/readthroughpombe/results/readthrough/peak_calling/2018_02_28_*_RT_forward.bed" \
  --bam "/scratch/lmodolo/readthroughpombe/results/readthrough/bams/metagene/*.bam" \
  -resume -w /scratch/lmodolo/work/ -with-dag /scratch/lmodolo/readthroughpombe/results/readthrough/metagene_dag.pdf \
  -with-timeline /scratch/lmodolo/readthroughpombe/results/readthrough/metagene_timeline

