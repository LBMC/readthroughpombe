#!/bin/bash

docker build src/func/docker -t 'pipeline:0.0.1' && \
bin/nextflow src/pipeline.nf -c src/pipeline.config -profile docker --fastq "data/examples/tiny_dataset/fastq/*tiny_S.fastq" --fasta "data/examples/tiny_dataset/fasta/tiny_v2.fasta" --todo "fastqc+cutadapt+urqt+fastqc+multiqc"
bin/nextflow src/pipeline.nf -c src/pipeline.config -profile docker --fastq "data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq" --fasta "data/examples/tiny_dataset/fasta/tiny_v2.fasta" --todo "fastqc+cutadapt+urqt+fastqc+multiqc"
