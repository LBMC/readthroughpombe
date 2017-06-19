#!/bin/sh
docker build src/pipe/quality_control -t 'quality_control:0.0.1'
bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'cutadapt' -c src/nextflow.config -with-docker quality_control:0.0.1
bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'cutadapt' -c src/nextflow.config -with-docker quality_control:0.0.1
bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'UrQt' -c src/nextflow.config -with-docker quality_control:0.0.1
bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'UrQt' -c src/nextflow.config -with-docker quality_control:0.0.1
