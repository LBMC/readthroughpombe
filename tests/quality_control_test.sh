#!/bin/sh
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'cutadapt' && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'cutadapt' && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'UrQt' && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'UrQt'
