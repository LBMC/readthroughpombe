#!/bin/sh
bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/*.fastq' --trimmer 'UrQt'
bin/nextflow src/pipe/quality_control.nf --fastq_files 'data/examples/*.fastq' --trimmer 'cutadapt'
