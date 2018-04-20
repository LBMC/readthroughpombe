# Methods

## data acquisition

Paragraph on the NGS sequencing of 4 strains

We downloaded the version 2.30 of *S. Pombe* genome in fasta format and the 
corresponding gff3 annotations on the ebi website (2017/03/22) (Supplementary 
section N).

## NGS analysis

All the NGS analysis scripts are available in the following git repository: giturl
We wrote most of these analyses as nextflow (v0.25.7) pipelines with Docker 
(v18.02.0-ce). In those pipelines, we used the following software’s:

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
the `src` folder of the git repository. The `init.sh` script contains 
instructions to get nextflow and initialize the dockers containers.

## quality control and strandness identification

The following analysis is described in the file `1_QC_mapping.sh`.
This first analysis uses the fastq files and the fasta file.

The first step of the analysis was to perform a quality control of the NGS files and to map them on the genome to determine their strand. The reads were first process with cutadapt (-a AGATCGGAAGAG -g CTCTTCCGATCT) to remove any 
remaining adaptors. Then the reads were processed with UrQt (--t 20). We ran
FastQC in command line mode before and after trimming and summarized its
output with MultiQC.

After the quality control, we built the reverse complement of the reads in 
the fastq files with seqkit.

Finally, we indexed the genome (-build) and mapped the trimmed fastq 
(--very-sensitive) with Bowtie2. We saved the mapping output in bam format with samtools (view -Sb).

## readthrough detection

The following analysis is described in the file `2_readthrough_detection.sh`.

We define a readthrough even as a genomic region where reads map in an intergenic
region after the 5’ end of a transcript annotation. In this section we used a
broad peak caller for chip-seq analysis to identify genomic regions where
enough reads are mapping to be detected as a peak. We then consider peak
annotation, starting after an annotated transcript, as putative readthrough 
annotation.

We ran the readthrough detection independently for the rrp6D and cut14-208
mutants. For each analysis the mutant was compared to the wild type.

Bam files were sorted and filtered as forward or reverse mapping with the
samtools tool (view -hb -F 0x10 or view -hb -f 0x10). The gff annotation file
was converted to bed with convert2bed from bedtools. Finally, the genome was
indexed with samtools (faidx) and the annotation file was split into a
forward and a reverse annotation with bedtools (complement).

To identify readthrough events, we started by removing any read
corresponding to a transcript with samtools view (-hb bams -L annotation). The
resulting bams were sorted with samtools sort. The bams and genome index were
given to Music to compute the genome mappability and preprocess the alignments.
The resulting sam files were sorted and dedupliated with Music.

We then ran Music (-get_per_win_p_vals_vs_FC -begin_l 50 -end_l 500 -step 1.1 
-l_mapp 50 -l_frag 363 -q_val 1 -l_p 0) to detect broad peaks, corresponding to
genomic regions were enough reads are found but which are not corresponding to
annotated transcripts. As Music gives a bed annotation file for each replicate,
we merged these annotations to only retain peaks found in the majority of
replicates (more than 2 for 3 replicates here). We then further filtered this
annotation by removing peaks whose starting position was more than 2 reads
length from the 5’ end of the nearest transcript on the same strand. We merged the transcripts and peaks annotation to associate a given peaks to the
nearest transcripts.

## readthrough quantification

The following analysis is described in the file `3_readthrough_quantification.sh`.

With the putative genomic region where readthrough event may be detected we
consider the problem of quantifying readthrough events as an alternative splicing
problem were mRNA fragments can correspond to the sort original form of the transcript or to long reathrough form of the transcript.

To have comparable analysis between rrp6D and cut14-208, we start by merging their respective readthrough annotations and sorting them with bedtools (sort).
The analysis was run separately for the forward and reverse strand of the genome.

We extracted the forward and reverse reads from the bam files obtained in the
readthrough detection analysis with bedtools (bamtofastq). We generated the
list of transcript sequences from the genome and the annotation with bedtools
(getfasta -s). The transcript sequences were then indexed with kallisto (index
-k 31 --make-unique). The quantification was done with kallisto (quant --single
-l 363.4 -s 85.53354). We ran the quantification separately on the transcript alone and on
the transcript with readthrough annotation.

## readthrough differential expression analysis
 
The following analysis is described in the file `4_readthrough_DEA.sh`.

With putative readthrough events annotated and quantified, we wanted to test if
the number of fragments corresponding to the readthrough events in the mutant was
significantly higher than in the wild type. For this we used the package DESeq2
(v1.16.1) with R (v3.4.4). We tested for a log2 fold change superior to 0,
compared to the wild type, to declare the read though form of a transcript 
present in the mutant. For the transcript deregulation analysis we tested for log2
fold change superior to 0.5 or inferior to -0.5. For all the analysis we selected
p-values with an FDR <= 0.05.





























