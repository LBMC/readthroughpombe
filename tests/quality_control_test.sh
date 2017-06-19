#!/bin/sh
# bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq.gz' --trimmer 'UrQt' --paired false -resume && \
# bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq.gz' --trimmer 'cutadapt' --paired false -resume && \
# bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq.gz' --trimmer 'UrQt' --paired true
docker build src/pipe/quality_control -t 'quality_control:0.0.1'
bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'UrQt' --paired true -c src/nextflow.config -with-docker quality_control:0.0.1
