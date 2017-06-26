#!/bin/sh
cp src/func/file_handle.py src/pipe/quality_control/ && \
docker build src/pipe/quality_control -t 'quality_control:0.0.1' && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'cutadapt' && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'cutadapt' && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'UrQt' && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_docker --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'UrQt'
