#!/bin/sh
# bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/*.fastq' --trimmer 'UrQt' --paired false -resume && \
# bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/*.fastq' --trimmer 'cutadapt' --paired false -resume && \
# bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/*_R{1,2}.fastq.gz' --trimmer 'UrQt' --paired true -resume
bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/*_R{1,2}.fastq.gz' --trimmer 'UrQt' --paired true -resume -c src/nextflow.config --with-docker 'src/pipe/mapping/'
