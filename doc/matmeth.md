# Methods

## data aquisition

Paragraph on the NGS sequencing of 4 strains

We downloaded the version 2.30 of *S. Pombe* genome in fasta format and the corresponding gff3 annotation on the ebi website (2017/03/22) (Supplementary section N).

## NGS analysis

All the NGS analysis scripts are available in the following git repository: giturl
We wrote most of these analysis as nextflow (v0.25.7) pipeline with Docker (v18.02.0-ce).
In those pipeline, we used the following softwares:

- FastQC (v0.11.5) and MultiQC (v0.9) for the quality control analysis,
- cutadapt (v1.14) to trim any remaining adaptators,
- UrQt (v1.0.17) to trim reads on their quality,
- Bowtie2 (v2.3.2) for the mapping,
- Kallisto (v0.43.1) for the quantification
- samtools (v1.5) and bedtools (v2.25.0) for various files transformation,
- Music (v6613c53) for the peak detection.

The bash scripts corresponding to the following subsections are available in the `src` folder of the git repository. The `init.sh` script contains the instructions to get nextflow and initialize the dockers containers.

## quality control and strandness identification

The following analysis is described in the file `1_QC_mapping.sh`.
This first analysis uses the fastq files and the fasta file.

The first step of the analysis was to perform a quality control of the NGS files and to map them on the genome to determine their strand. The reads were first process with cutadapt (-a AGATCGGAAGAG -g CTCTTCCGATCT) to remove any remaining adaptors. Then the reads were processed with UrQt (--t 20). We ran FastQC in commande line mode before and after trimming and summarized it's output with MultiQC.
After the quality control, we indexed the genome (-build) and mapped the trimmed fastq (--very-sensitive) with bowtie2.We saved the mapping output in bam format with samtools (view -Sb).

## readthrough detection

The following analysis is described in the file `2_readthrough_detection.sh`.



## readthrough analysis
