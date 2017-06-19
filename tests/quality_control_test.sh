#!/bin/sh
bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'cutadapt' -c src/nextflow.config
bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'cutadapt' -c src/nextflow.config
bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq' --trimmer 'UrQt' -c src/nextflow.config
bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/tiny_dataset/fastq/*tiny_S.fastq' --trimmer 'UrQt' -c src/nextflow.config
