#!/bin/bash

docker build src/func/docker -t 'pipeline:0.0.1' && \
bin/nextflow src/pipeline.nf -c src/pipeline.config -profile docker --fastq "data/examples/tiny_dataset/fastq/*tiny_S.fastq" --fasta "data/examples/tiny_dataset/fasta/tiny_v2.fasta" --annot "data/examples/tiny_dataset/annot/tiny.gff" --todo "fastqc+cutadapt+urqt+fastqc+multiqc+split_ref"
bin/nextflow src/pipeline.nf -c src/pipeline.config -profile docker --fastq "data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq" --fasta "data/examples/tiny_dataset/fasta/tiny_v2.fasta" --annot "data/examples/tiny_dataset/annot/tiny.gff" --todo "fastqc+cutadapt+urqt+fastqc+multiqc+split_ref"
