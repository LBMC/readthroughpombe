bin/nextflow src/RNASeq/src/pipe/quality_control.nf -c src/RNASeq/src/nextflow.config -profile quality_control --fastq_files "data/RNASeq/*fastq.gz" --cpu 2 --fastqc "/usr/local/bin/FastQC/fastqc" --trimmer "cutadapt" --urqt "/usr/local/bin/UrQt/UrQt" --paired false --quality_threshold 20 --pigz_version "2.3" --do_adapter_removal false

