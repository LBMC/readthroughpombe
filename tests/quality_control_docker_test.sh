#!/bin/sh
cp src/file_handle/src/file_handle.py src/pipe/quality_control/ && \
docker build src/pipe/quality_control -t 'quality_control:0.0.1' && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'cutadapt' --cpu 2 && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'cutadapt' --cpu 2  && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'UrQt' --cpu 2  && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'UrQt' --cpu 2 
