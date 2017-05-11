#!/bin/sh
bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/*.fastq' --trimmer 'UrQt' && \
ls -lR results/quality_control/ && \
bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/*.fastq' --trimmer 'cutadapt' && \
ls -lR results/quality_control/
