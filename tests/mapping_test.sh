#!/bin/sh
# bin/nextflow src/pipe/mapping.nf --fastq_files 'data/examples/*tiny_R{1,2}.fastq.gz' --mapper 'salmon' --paired true --reference 'data/examples/*tiny.fasta.gz' --annotation 'data/examples/*tiny.gff'
# bin/nextflow src/pipe/mapping.nf --fastq_files 'data/examples/*tiny_S.fastq.gz' --mapper 'salmon' --paired false --reference 'data/examples/*tiny.fasta.gz' --annotation 'data/examples/*tiny.gff'
# bin/nextflow src/pipe/mapping.nf --fastq_files 'data/examples/*tiny_R{1,2}.fastq.gz' --mapper 'kallisto' --paired true --reference 'data/examples/*tiny.fasta.gz' --annotation 'data/examples/*tiny.gff'
# bin/nextflow src/pipe/mapping.nf --fastq_files 'data/examples/*tiny_S.fastq.gz' --mapper 'kallisto' --paired false --reference 'data/examples/*tiny.fasta.gz' --annotation 'data/examples/*tiny.gff'
# bin/nextflow src/pipe/mapping.nf --fastq_files 'data/examples/*tiny_R{1,2}.fastq.gz' --mapper 'bowtie2' --paired true --reference 'data/examples/*tiny.fasta.gz' --annotation 'data/examples/*tiny.gff'
# bin/nextflow src/pipe/mapping.nf --fastq_files 'data/examples/*tiny_S.fastq.gz' --mapper 'bowtie2' --paired false --reference 'data/examples/*tiny.fasta.gz' --annotation 'data/examples/*tiny.gff'
bin/nextflow src/pipe/mapping.nf --fastq_files 'data/examples/*tiny_R{1,2}.fastq.gz' --mapper 'bowtie2' --paired true --reference 'data/examples/*tiny.fasta.gz' --annotation 'data/examples/*tiny.gff' --quantifier 'rsem'
