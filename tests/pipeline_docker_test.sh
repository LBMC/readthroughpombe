#!/bin/bash

docker build src/func/mapping -t 'mapping:0.0.1'&& \
docker build src/func/mapping -t 'quality_control:0.0.1' && \
bin/nextflow src/pipeline.nf -c src/pipeline.config -profile docker --fastq "data/examples/tiny_dataset/fastq/*tiny_S.fastq" --todo "fastqc+cutadapt+urqt" -with-docker quality_control:0.0.1
