#!/usr/bin/env nextflow

/*
* Copyright laurent modolo for the LBMC UMR 5239 ©.
* contributor(s) : laurent modolo (2017)
*
* laurent.modolo@ens-lyon.fr
*
* This software is a computer program whose purpose is to detect
* readthrough events from stranded RNASeq mapped reads
*
* This software is governed by the CeCILL  license under French law and
* abiding by the rules of distribution of free software.  You can  use,
* modify and/ or redistribute the software under the terms of the CeCILL
* license as circulated by CEA, CNRS and INRIA at the following URL
* "http://www.cecill.info".
*
* As a counterpart to the access to the source code and  rights to copy,
* modify and redistribute granted by the license, users are provided only
* with a limited warranty  and the software's author,  the holder of the
* economic rights,  and the successive licensors  have only  limited
* liability.
*
* In this respect, the user's attention is drawn to the risks associated
* with loading,  using,  modifying and/or developing or reproducing the
* software by the user in light of its specific status of free software,
* that may mean  that it is complicated to manipulate,  and  that  also
* therefore means  that it is reserved for developers  and  experienced
* professionals having in-depth computer knowledge. Users are therefore
* encouraged to load and test the software's suitability as regards their
* requirements in conditions enabling the security of their systems and/or
* data to be ensured and,  more generally, to use and operate it in the
* same conditions as regards security.
*
* The fact that you are presently reading this means that you have had
* knowledge of the CeCILL license and that you accept its terms.
*/

params.verbose = false
params.annotation_forward = ""
params.annotation_reverse = ""
params.fasta = ""
params.bam = ""
params.mean_size = 200
params.sd_size = 200

bam_files = Channel
  .fromPath(params.bam)
  .ifEmpty {
    exit 1, "Cannot find any bam file matching: ${params.bam}"
  }

genome_file = Channel
  .fromPath(params.fasta)
  .ifEmpty {
    exit 1, "Cannot find any fasta file matching: ${params.fasta}"
  }

annotation_forward_files = Channel
  .fromPath(params.annotation_forward)
  .ifEmpty {
    exit 1, "Cannot find any annotation file matching: ${params.annotation_forward}"
  }

annotation_reverse_files = Channel
  .fromPath(params.annotation_reverse)
  .ifEmpty {
    exit 1, "Cannot find any annotation file matching: ${params.annotation_forward}"
  }

process merge_annotation_forward {
  echo params.verbose
  publishDir "results/readthrough/transcript/", mode: 'copy'
  cpus 1
  input:
    file annotation from annotation_forward_files.collect()
  output:
    file "*annotation_forward.bed" into rt_forward_bed
    file "*transcript_forward.bed" into t_forward_bed
  script:
    """
    cat ${annotation} > annotation_merge.bed

    grep "ID=transcript:" annotation_merge.bed | \
      grep -v "transcript_RT" | \
      awk '{if(\$2 < \$3) {print \$0}}' > \
      annotation_merge_T.bed

    grep "transcript_RT" annotation_merge.bed | \
      grep "ID=transcript:" | \
      awk '{if(\$2 < \$3) {print \$0}}' | \
      sed 's/\\(.*;Name=\\)\\(.*\\)\\(;Parent.*\\)/\\1\\2_RT\\3/g' | \
      sed 's/transcript:/transcript:RT_/g' > \
      annotation_merge_RT.bed

    cat annotation_merge_T.bed annotation_merge_RT.bed | \
      sort -k4 | uniq | \
      bedtools sort -i stdin > annotation_forward.bed

    bedtools sort -i annotation_merge_T > transcript_forward.bed

    file_handle.py -f annotation_forward.bed transcript_forward.bed
    """
}

process merge_annotation_reverse {
  echo params.verbose
  publishDir "results/readthrough/transcript/", mode: 'copy'
  cpus 1
  input:
    file annotation from annotation_reverse_files.collect()
  output:
    file "*annotation_reverse.bed" into rt_reverse_bed
    file "*transcript_reverse.bed" into t_reverse_bed
  script:
    """
    cat ${annotation} > annotation_merge.bed

    grep "ID=transcript:" annotation_merge.bed | \
      grep -v "transcript_RT" | \
      awk '{if(\$2 < \$3) {print \$0}}' > \
      annotation_merge_T.bed

    grep "transcript_RT" annotation_merge.bed | \
      grep "ID=transcript:" | \
      awk '{if(\$2 < \$3) {print \$0}}' | \
      sed 's/\\(.*;Name=\\)\\(.*\\)\\(;Parent.*\\)/\\1\\2_RT\\3/g' | \
      sed 's/transcript:/transcript:RT_/g' > \
      annotation_merge_RT.bed

    cat annotation_merge_T.bed annotation_merge_RT.bed | \
      sort -k4 | uniq | \
      bedtools sort -i stdin > annotation_reverse.bed

      bedtools sort -i annotation_merge_T > transcript_reverse.bed

      file_handle.py -f annotation_reverse.bed transcript_reverse.bed
    """
}

process sort_bam {
  echo params.verbose
  cpus 4
  input:
    file bams from bam_files
  output:
    file "*_sorted.bam" into sorted_bams
  script:
    """
    samtools sort -@ ${task.cpus} -O BAM -o ${bams}_sorted.bam ${bams}
    find . -name "*_sorted.bam" | sed 's/\\.bam_sorted\\.bam//g' | \
    awk '{system("mv "\$0".bam_sorted.bam "\$0"_sorted.bam")}'
    find . -name "*_sorted.bam" |\
    sed 's/\\(.*\\)_rev.*/\\1/g' |\
    sed 's/\\(.*\\)_trim.*/\\1/g' |\
    awk '{system("mv "\$0"*_sorted.bam "\$0"_sorted.bam")}'
    file_handle.py -f *.bam
    """
}

process split_bam {
  echo params.verbose
  cpus 2
  input:
    file bams from sorted_bams
  output:
    file "forward/*_forward.bam*" into forward_bams
    file "reverse/*_reverse.bam*" into reverse_bams
  script:
  """
    samtools view -hb -F 0x10 ${bams} > ${bams}_forward.bam &
    samtools view -hb -f 0x10 ${bams} > ${bams}_reverse.bam
    mkdir -p forward reverse
    find . -name "*_forward.bam" | sed 's/\\.bam_forward\\.bam//g' | \
    awk '{system("mv "\$0".bam_forward.bam "\$0"_forward.bam")}'
    find . -name "*_reverse.bam" | sed 's/\\.bam_reverse\\.bam//g' | \
    awk '{system("mv "\$0".bam_reverse.bam "\$0"_reverse.bam")}'
    file_handle.py -f *_forward.bam *_reverse.bam
    mv *_forward.bam forward/
    mv *_reverse.bam reverse/
  """
}

process bam_2_fastq_forward {
  echo params.verbose
  publishDir "results/readthrough/fastq/forward/", mode: 'copy'
  cpus 1
  input:
    file bams from forward_bams
  output:
    file "*_forward.fastq" into fastq_forward
  script:
  """
    bedtools bamtofastq -i ${bams} -fq ${bams}.fastq
    find . -name "*.bam.fastq" | sed 's/\\.bam\\.fastq//g' | \
    awk '{system("mv "\$0".bam.fastq "\$0".fastq")}'
    file_handle.py -f *_forward.fastq
  """
}

process bam_2_fastq_reverse {
  echo params.verbose
  publishDir "results/readthrough/fastq/reverse/", mode: 'copy'
  cpus 1
  input:
    file bams from reverse_bams
  output:
    file "*_reverse.fastq" into fastq_reverse
  script:
  """
    bedtools bamtofastq -i ${bams} -fq ${bams}.fastq
    find . -name "*.bam.fastq" | sed 's/\\.bam\\.fastq//g' | \
    awk '{system("mv "\$0".bam.fastq "\$0".fastq")}'
    file_handle.py -f *_reverse.fastq
  """
}

genome_file.into{
  genome_file_rt_forward;
  genome_file_rt_reverse;
  genome_file_t_forward;
  genome_file_t_reverse
}

process rt_forward {
  echo params.verbose
  cpus 1
  input:
    file genome from genome_file_rt_forward
    file annotation from rt_forward_bed
  output:
    file "*.fasta" into rt_forward_fasta
  script:
  """
    bedtools getfasta -s -name -fi ${genome} -bed ${annotation} -fo rt_forward.fasta
    file_handle.py -f rt_forward.fasta
  """
}

process rt_reverse {
  echo params.verbose
  publishDir "results/readthrough/transcript/", mode: 'copy'
  cpus 1
  input:
    file genome from genome_file_rt_reverse
    file annotation from rt_reverse_bed
  output:
    file "*.fasta" into rt_reverse_fasta
  script:
  """
    bedtools getfasta -s -name -fi ${genome} -bed ${annotation} -fo rt_reverse.fasta
    file_handle.py -f rt_reverse.fasta
  """
}

process t_forward {
  echo params.verbose
  cpus 1
  input:
    file genome from genome_file_t_forward
    file annotation from t_forward_bed
  output:
    file "*.fasta" into t_forward_fasta
  script:
  """
    bedtools getfasta -s -name -fi ${genome} -bed ${annotation} -fo t_forward.fasta
    file_handle.py -f t_forward.fasta
  """
}

process t_reverse {
  echo params.verbose
  publishDir "results/readthrough/transcript/", mode: 'copy'
  cpus 1
  input:
    file genome from genome_file_t_reverse
    file annotation from t_reverse_bed
  output:
    file "*.fasta" into t_reverse_fasta
  script:
  """
    bedtools getfasta -s -name -fi ${genome} -bed ${annotation} -fo t_reverse.fasta
    file_handle.py -f t_reverse.fasta
  """
}

process rt_indexing_forward {
  echo params.verbose
  publishDir "results/readthrough/transcript/", mode: 'copy'
  cpus 1
  input:
    file genome from rt_forward_fasta
  output:
    file "*index*" into rt_forward_index
  script:
  """
    kallisto index -k 31 --make-unique -i forward.index ${genome} &> kallisto_rt_index_forward_report.txt
    file_handle.py -f *index*
  """
}

process rt_indexing_reverse {
  echo params.verbose
  cpus 1
  input:
    file genome from rt_reverse_fasta
  output:
    file "*index*" into rt_reverse_index
  script:
  """
    kallisto index -k 31 --make-unique -i reverse.index ${genome} &> kallisto_rt_index_reverse_report.txt
    file_handle.py -f *index*
  """
}

process t_indexing_forward {
  echo params.verbose
  publishDir "results/readthrough/transcript/", mode: 'copy'
  cpus 1
  input:
    file genome from t_forward_fasta
  output:
    file "*index*" into t_forward_index
  script:
  """
    kallisto index -k 31 --make-unique -i forward.index ${genome} &> kallisto_t_index_forward_report.txt
    file_handle.py -f *index*
  """
}

process t_indexing_reverse {
  echo params.verbose
  cpus 1
  input:
    file genome from t_reverse_fasta
  output:
    file "*index*" into t_reverse_index
  script:
  """
    kallisto index -k 31 --make-unique -i reverse.index ${genome} &> kallisto_t_index_reverse_report.txt
    file_handle.py -f *index*
  """
}

fastq_forward.into{
  fastq_forward_rt;
  fastq_forward_t
}

fastq_reverse.into{
  fastq_reverse_rt;
  fastq_reverse_t
}

process quantification_rt_forward {
  echo params.verbose
  cpus 4
  publishDir "results/readthrough/quantification/forward/", mode: 'copy'
  input:
    file index from rt_forward_index.collect()
    file fastq from fastq_forward_rt
  output:
    file "*" into quantification_forward_rt
  script:
  tagname = (fastq =~ /(.*\/){0,1}d{0,1}(.*)\.fast[aq](\.gz){0,1}/)[0][2]
  """
    kallisto quant -i *forward.index -t ${task.cpus} --single -l ${params.mean_size} \
      -s ${params.sd_size} -o ./ ${fastq} &> ${tagname}_rt_kallisto_report.txt
    mv abundance.tsv ${tagname}_rt.tsv
    mv run_info.json ${tagname}_info_rt.json
    mv abundance.h5 ${tagname}_rt.h5
    file_handle.py -f *
  """
}

process quantification_rt_reverse {
  echo params.verbose
  cpus 4
  publishDir "results/readthrough/quantification/reverse/", mode: 'copy'
  input:
    file index from rt_reverse_index.collect()
    file fastq from fastq_reverse_rt
  output:
    file "*" into quantification_reverse_rt
  script:
  tagname = (fastq =~ /(.*\/){0,1}d{0,1}(.*)\.fast[aq](\.gz){0,1}/)[0][2]
  """
    kallisto quant -i *reverse.index -t ${task.cpus} --single -l ${params.mean_size} \
      -s ${params.sd_size} -o ./ ${fastq} &> ${tagname}_rt_kallisto_report.txt
    mv abundance.tsv ${tagname}_rt.tsv
    mv run_info.json ${tagname}_info_rt.json
    mv abundance.h5 ${tagname}_rt.h5
    file_handle.py -f *
  """
}

process quantification_t_forward {
  echo params.verbose
  cpus 4
  publishDir "results/readthrough/quantification/forward/", mode: 'copy'
  input:
    file index from t_forward_index.collect()
    file fastq from fastq_forward_t
  output:
    file "*" into quantification_forward_t
  script:
  tagname = (fastq =~ /(.*\/){0,1}d{0,1}(.*)\.fast[aq](\.gz){0,1}/)[0][2]
  """
    kallisto quant -i *forward.index -t ${task.cpus} --single -l ${params.mean_size} \
      -s ${params.sd_size} -o ./ ${fastq} &> ${tagname}_t_kallisto_report.txt
    mv abundance.tsv ${tagname}_t.tsv
    mv run_info.json ${tagname}_info_t.json
    mv abundance.h5 ${tagname}_t.h5
    file_handle.py -f *
  """
}

process quantification_t_reverse {
  echo params.verbose
  cpus 4
  publishDir "results/readthrough/quantification/reverse/", mode: 'copy'
  input:
    file index from t_reverse_index.collect()
    file fastq from fastq_reverse_t
  output:
    file "*" into quantification_reverse_t
  script:
  tagname = (fastq =~ /(.*\/){0,1}d{0,1}(.*)\.fast[aq](\.gz){0,1}/)[0][2]
  """
    kallisto quant -i *reverse.index -t ${task.cpus} --single -l ${params.mean_size} \
      -s ${params.sd_size} -o ./ ${fastq} &> ${tagname}_t_kallisto_report.txt
    mv abundance.tsv ${tagname}_t.tsv
    mv run_info.json ${tagname}_info_t.json
    mv abundance.h5 ${tagname}_t.h5
    file_handle.py -f *
  """
}
