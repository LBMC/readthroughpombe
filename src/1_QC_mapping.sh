# script to run the pipeline.nf on the RNASeq data

# with docker

# QC analysis
bin/nextflow src/pipeline.nf -c src/pipeline.config -profile docker --fastq "data/RNASeq/trimming/2017_04_14*.fastq.gz" --todo "fastq+cutadapt+urqt+fastqc+multiqc"

# mapping
bin/nextflow src/pipeline.nf -c src/pipeline.config -profile docker --fastq "results/quality_control/trimming/2017_09_*.fastq.gz" --fasta 'data/2017_09_19_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fasta.gz' --todo "bowtie2"
