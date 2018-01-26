# script to run the pipeline readthrough.nf on the data

bin/nextflow src/readthrough.nf -c src/pipe/conf/readthrough_docker.config --name "rrp6D" --bams "data/bams_rrp6D.csv" --annotation "data/2017_09_19_schizosaccharomyces_pombe.chr.gff3" --index "results/mapping/indexing/2018_01_19_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.index*" --genome "data/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa" --reads_size 50 --frag_size 363

bin/nextflow src/readthrough.nf -c src/pipe/conf/readthrough_docker.config --name "cut14" --bams "data/bams_cut14.csv" --annotation "data/2017_09_19_schizosaccharomyces_pombe.chr.gff3" --index "results/mapping/indexing/2018_01_19_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.index*" --genome "data/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa" --reads_size 50 --frag_size 363 -resume
