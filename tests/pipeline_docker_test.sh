#!/bin/bash

docker build src/func/docker -t 'pipeline:0.0.1' && \
bin/nextflow src/pipeline.nf -c src/pipeline.config -profile docker --fastq "data/examples/tiny_dataset/fastq/*tiny_S.fastq" --todo "fastqc+cutadapt+urqt" -with-docker pipeline:0.0.1
