#!/bin/sh
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'cutadapt' --cpu 2 && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'cutadapt' --cpu 2 && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'UrQt' --cpu 2 && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'UrQt' --cpu 2

ls -l data/examples/tiny_dataset/fastq/*.fastq | awk '{system("gzip "$9" -c > "$9".gz")}'

bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq.gz' --trimmer 'cutadapt' --cpu 2 && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq.gz' --trimmer 'cutadapt' --cpu 2 && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq.gz' --trimmer 'UrQt' --cpu 2 && \
bin/nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq.gz' --trimmer 'UrQt' --cpu 2
