#!/bin/sh
bin/nextflow src/pipe/mapping.nf --fastq_files 'data/examples/*reads_{1,2}.fastq.gz' --mapper 'salmon' --paired true --reference 'data/examples/*transcripts.fasta.gz'
bin/nextflow src/pipe/mapping.nf --fastq_files 'data/examples/*reads_1.fastq.gz' --mapper 'salmon' --paired false --reference 'data/examples/*transcripts.fasta.gz'
bin/nextflow src/pipe/mapping.nf --fastq_files 'data/examples/*reads_{1,2}.fastq.gz' --mapper 'kallisto' --paired true --reference 'data/examples/*transcripts.fasta.gz'
bin/nextflow src/pipe/mapping.nf --fastq_files 'data/examples/*reads_1.fastq.gz' --mapper 'kallisto' --paired false --reference 'data/examples/*transcripts.fasta.gz'
