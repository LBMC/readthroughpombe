# script to run the pipeline readthrough.nf on the data

bin/nextflow src/readthrough.nf -c src/pipe/conf/readthrough_docker.config -profile docker --bams "results/mapping/mapping/2018_04_14*.bam" --annotation "data/2017_09_19_schizosaccharomyces_pombe.chr.gff3" -resume
