#!/bin/bash

cp data/examples/tiny_dataset/fastq/tiny_R1.fastq data/examples/tiny_dataset/fastq/c1tiny_R1.fastq && \
cp data/examples/tiny_dataset/fastq/tiny_R2.fastq data/examples/tiny_dataset/fastq/c1tiny_R2.fastq && \
cp data/examples/tiny_dataset/fastq/tiny_R1.fastq data/examples/tiny_dataset/fastq/c2tiny_R1.fastq && \
cp data/examples/tiny_dataset/fastq/tiny_R2.fastq data/examples/tiny_dataset/fastq/c2tiny_R2.fastq && \
perl -pi -e 's/queue.*/queue = "E5_test"/g' src/pipe/conf/pipeline_sge.config
perl -pi -e 's/cpus.*/cpus = 1/g' src/pipe/conf/pipeline_sge.config
perl -pi -e 's/time.*//g' src/pipe/conf/pipeline_sge.config
perl -pi -e 's/penv.*//g' src/pipe/conf/pipeline_sge.config


module load nextflow/0.25.1

nextflow src/pipeline.nf -c src/pipeline.config -profile sge --fastq "data/examples/tiny_dataset/fastq/*tiny_S.fastq" --fasta "data/examples/tiny_dataset/fasta/tiny_v2.fasta" --annot "data/examples/tiny_dataset/annot/tiny.gff" --todo "fastqc+cutadapt+cutadapt+fastqc+multiqc+split_ref+kallisto"
nextflow src/pipeline.nf -c src/pipeline.config -profile sge --fastq "data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq" --fasta "data/examples/tiny_dataset/fasta/tiny_v2.fasta" --annot "data/examples/tiny_dataset/annot/tiny.gff" --todo "fastqc+cutadapt+cutadapt+fastqc+multiqc+split_ref+kallisto"

nextflow src/pipeline.nf -c src/pipeline.config -profile sge --fastq "data/examples/tiny_dataset/fastq/*tiny_S.fastq" --fasta "data/examples/tiny_dataset/fasta/tiny_v2.fasta" --annot "data/examples/tiny_dataset/annot/tiny.gff" --todo "fastqc+cutadapt+cutadapt+fastqc+multiqc+bowtie2+htseq"
nextflow src/pipeline.nf -c src/pipeline.config -profile sge --fastq "data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq" --fasta "data/examples/tiny_dataset/fasta/tiny_v2.fasta" --annot "data/examples/tiny_dataset/annot/tiny.gff" --todo "fastqc+cutadapt+cutadapt+fastqc+multiqc+split_ref+bowtie2+htseq"

nextflow src/pipeline.nf -c src/pipeline.config -profile sge --fastq "data/examples/tiny_dataset/fastq/*tiny_S.fastq" --fasta "data/examples/tiny_dataset/fasta/tiny_v2.fasta" --annot "data/examples/tiny_dataset/annot/tiny.gff" --todo "fastqc+cutadapt+cutadapt+fastqc+multiqc+bowtie2+rsem"
nextflow src/pipeline.nf -c src/pipeline.config -profile sge --fastq "data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq" --fasta "data/examples/tiny_dataset/fasta/tiny_v2.fasta" --annot "data/examples/tiny_dataset/annot/tiny.gff" --todo "fastqc+cutadapt+cutadapt+fastqc+multiqc+bowtie2+rsem"

nextflow src/pipeline.nf -c src/pipeline.config -profile sge --fastq "data/examples/tiny_dataset/fastq/*tiny_S.fastq" --fasta "data/examples/tiny_dataset/fasta/tiny_v2.fasta" --annot "data/examples/tiny_dataset/annot/tiny.gff" --todo "fastqc+cutadapt+cutadapt+fastqc+multiqc+kallisto"
nextflow src/pipeline.nf -c src/pipeline.config -profile sge --fastq "data/examples/tiny_dataset/fastq/*tiny_R{1,2}.fastq" --fasta "data/examples/tiny_dataset/fasta/tiny_v2.fasta" --annot "data/examples/tiny_dataset/annot/tiny.gff" --todo "fastqc+cutadapt+cutadapt+fastqc+multiqc+kallisto"

git checkout src/pipe/conf/pipeline_sge.config
