# Methods

## data acquisition

Paragraph on the NGS sequencing of 4 strains

We downloaded the version 2.30 of *S. Pombe* genome in fasta format and the 
corresponding gff3 annotation on the ebi website (2017/03/22) (Supplementary 
section N).

## NGS analysis

All the NGS analysis scripts are available in the following git repository: giturl
We wrote most of these analysis as nextflow (v0.25.7) pipelines with Docker 
(v18.02.0-ce). In those pipeline, we used the following software's:

- FastQC (v0.11.5) and MultiQC (v0.9) for the quality control analysis,
- cutadapt (v1.14) to trim any remaining adaptors,
- UrQt (v1.0.17) to trim reads on their quality,
- Bowtie2 (v2.3.2) for the mapping,
- Kallisto (v0.43.1) for the quantification
- samtools (v1.5) and bedtools (v2.25.0) for various files transformation,
- Music (v6613c53) for the peak detection.
- segkit (v0.7.2) to compute reverse complement of sequences
- bedtools (2.25.0) for annotation and alignment processing

The bash scripts corresponding to the following subsections are available in 
the `src` folder of the git repository. The `init.sh` script contains the 
instructions to get nextflow and initialize the dockers containers.

## quality control and strandness identification

The following analysis is described in the file `1_QC_mapping.sh`.
This first analysis uses the fastq files and the fasta file.

The first step of the analysis was to perform a quality control of the NGS 
files and to map them on the genome to determine their strand. The reads were 
first process with cutadapt (-a AGATCGGAAGAG -g CTCTTCCGATCT) to remove any 
remaining adaptors. Then the reads were processed with UrQt (--t 20). We ran 
FastQC in command line mode before and after trimming and summarized it's 
output with MultiQC.

After the quality control, we builed the reverse complement of the reads in 
the fastq files with seqkit.

Finaly, we indexed the genome (-build) and mapped the trimmed fastq 
(--very-sensitive) with Bowtie2. We saved the mapping output in bam format 
with samtools (view -Sb).

## readthrough detection

The following analysis is described in the file `2_readthrough_detection.sh`.

We can the readthrough detection independently for the rrp6D and cut14-208
mutants. For each analysis the mutant was compared to the wild type.

The bam files were sorted and filtered as forward or reverse mapping with the
samtools tool (view -hb -F 0x10 or view -hb -f 0x10). The gff annotation file
was comverted to bed with convert2bed from bedtools. Finaly, the genome was
indexed with samtools (faidx) to split the annotation file into a forward and a
reverse annotation with bedtools (complement).

To identify readthrough events were the reads associated with a genes continue
after the 5' end of his annotation, we started by removing any read
corresponding to a transcript with samtools view (-hb bams -L annotation). The
resulting bams were sorted with samtools sort. The bams and genome index were
givent to Music to compute the genome mappability and preprocess the alignments.
The resultsing sam file were sorted and dedupliated with Music.

We then ran Music (-get_per_win_p_vals_vs_FC -begin_l 50 -end_l 500 -step 1.1 
-l_mapp 50 -l_frag 363 -q_val 1 -l_p 0) to detect broad peaks, corresponding to
genomic regions were enought reads are found but which are not corresponding to
annotated transcripts. As Music gives a bed annotation file for each replicate,
we merged these annotation to only retain peaks found in the majority of
replicates (more than 2 for 3 replicates here). We then further filtered this
annotation by removing peaks whose starting position was at more than 2 reads
size from the 5' of the nearest transcript in the same reading frame. We merged
the transcripts and peaks annotation to associate a given peaks to the nearest
transcripts.

## readthrough analysis
