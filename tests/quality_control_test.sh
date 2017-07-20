#!/bin/sh
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'cutadapt' --cpu 2 && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'cutadapt' --cpu 2 && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'UrQt' --cpu 2 && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'UrQt' --cpu 2 
