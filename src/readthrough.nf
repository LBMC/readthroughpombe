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
params.bams = ""
params.annotation = ""
params.index = ""
params.genome = ""
params.reads_size = 0
params.frag_size = 200

bam_files = Channel.fromPath(params.bams)
  .ifEmpty { exit 1, "Cannot find any bams files matching: ${params.bams}" }
annotation_file = Channel.fromPath(params.annotation)
  .ifEmpty {
    exit 1, "Cannot find any annotation file matching: ${params.annotation}"
  }
index_file = Channel.fromPath(params.index)
  .ifEmpty {
    exit 1, "Cannot find any bowtie2 index_file file matching: ${params.index}"
  }
genome_file = Channel.fromPath(params.genome)
  .ifEmpty {
    exit 1, "Cannot find any fasta file matching: ${params.genome}"
  }

process sort_bam {
  echo params.verbose
  publishDir "results/readthrough/bam_sorting", mode: 'copy'
  cpus = 4
  input:
    file bam from bam_files
  output:
    file "*_sorted.bam" into sorted_bams
  script:
    """
    samtools sort -@ ${task.cpus} -O BAM -o ${bam}_sorted.bam ${bam}
    find . -name "*_sorted.bam" | sed 's/\\.bam_sorted\\.bam//g' | \
    awk '{system("mv "\$0".bam_sorted.bam "\$0"_sorted.bam")}'
    file_handle.py -f *.bam*
    """
}

process split_bam {
  echo params.verbose
  publishDir "results/readthrough/bam_spliting", mode: 'copy'
  cpus = 2
  input:
    file bam from sorted_bams
  output:
    file "*_forward.bam*" into forward_bams
    file "*_reverse.bam*" into reverse_bams
  script:
  """
    samtools view -hb -F 0x10 ${bam} > ${bam}_forward.bam &
    samtools view -hb -f 0x10 ${bam} > ${bam}_reverse.bam
    find . -name "*_forward.bam" | sed 's/\\.bam_forward\\.bam//g' | \
    awk '{system("mv "\$0".bam_forward.bam "\$0"_forward.bam")}'
    find . -name "*_reverse.bam" | sed 's/\\.bam_reverse\\.bam//g' | \
    awk '{system("mv "\$0".bam_reverse.bam "\$0"_reverse.bam")}'
    file_handle.py -f *_forward.bam *_reverse.bam
  """
}

process gff_to_bed {
  echo params.verbose
  publishDir "results/readthrough/annotation", mode: 'copy'
  cpus 4
  input:
    file annotation from annotation_file
  output:
    file "*_forward.bed*" into annotation_forward
    file "*_reverse.bed*" into annotation_reverse
  script:
  """
  ls -l
  convert2bed --input=gff --output=bed < ${annotation} > ${annotation}.bed
  awk '{if(\$6 == "+"){print \$0 >> "${annotation}_forward.bed"}; if(\$6 == "-"){print \$0 >> "${annotation}_reverse.bed"}'} ${annotation}.bed
  find . -name "*gff3_forward.bed" | sed 's/\\.gff3_forward\\.bed//g' | \
  awk '{system("mv "\$0".gff3_forward.bed "\$0"_forward.bed")}'
  find . -name "*_reverse.bed" | sed 's/\\.gff3_reverse\\.bed//g' | \
  awk '{system("mv "\$0".gff3_reverse.bed "\$0"_reverse.bed")}'
  ls -l
  """
}

process filter_forward {
  echo params.verbose
  publishDir "results/readthrough/bam_filtering", mode: 'copy'
  cpus 4
  input:
    file bam from forward_bams
    file annotation from annotation_forward.first()
  output:
    file "*_noannot.bam*" into filtered_forward_bams
  script:
  """
    samtools view -@ ${task.cpus} -hb ${bam} -L ${annotation} -U ${bam}_noannot.bam > /dev/null
    samtools index ${bam}_noannot.bam
    find . -name "*_noannot.bam" | sed 's/\\.bam_noannot\\.bam//g' | \
    awk '{system("mv "\$0".bam_noannot.bam "\$0"_noannot.bam")}'
    find . -name "*_noannot.bam.bai" | sed 's/\\.bam_noannot\\.bam.bai//g' | \
    awk '{system("mv "\$0".bam_noannot.bam.bai "\$0"_noannot.bam.bai")}'
    file_handle.py -f *_noannot.bam*
  """
}

process filter_reverse {
  echo params.verbose
  publishDir "results/readthrough/bam_filtering", mode: 'copy'
  cpus 4
  input:
    file bam from reverse_bams
    file annotation from annotation_reverse.first()
  output:
    file "*_noannot.bam*" into filtered_reverse_bams
  script:
  """
    samtools view -@ ${task.cpus} -hb ${bam} -L ${annotation} -U ${bam}_noannot.bam > /dev/null
    samtools index ${bam}_noannot.bam
    find . -name "*_noannot.bam" | sed 's/\\.bam_noannot\\.bam//g' | \
    awk '{system("mv "\$0".bam_noannot.bam "\$0"_noannot.bam")}'
    find . -name "*_noannot.bam.bai" | sed 's/\\.bam_noannot\\.bam.bai//g' | \
    awk '{system("mv "\$0".bam_noannot.bam.bai "\$0"_noannot.bam.bai")}'
    file_handle.py -f *_noannot.bam*
  """
}

process compute_mappability {
  echo params.verbose
  publishDir "results/readthrough/mappability", mode: copy
  cpu 4
  input:
    file index from index_files.collect()
    file genome from genome_file
  output:
    file "*.bin" into mappability
  script:
  """
  generate_multimappability_signal.csh ${genome} ${params.reads_size} ./
  bash temp_map_reads.csh
  bash temp_process_mapping.csh
  """
}

process music_preprocess_forward {
  echo params.verbose
  publishDir "results/readthrough/preprocessed", mode: copy
  cpu 1
  input:
    file bam from filtered_forward_bams
  output:
    file "preprocessed/*" into preprocessed_forward_bams
  script:
  """
  mkdir preprocessed
  samtools view ${bam} | \
  MUSIC -preprocess SAM stdin preprocessed/
  """
}

process music_preprocess_reverse {
  echo params.verbose
  publishDir "results/readthrough/preprocessed", mode: copy
  cpu 1
  input:
    file bam from filtered_reverse_bams
  output:
    file "preprocessed/*" into preprocessed_reverse_bams
  script:
  """
  mkdir preprocessed
  samtools view ${bam} | \
  MUSIC -preprocess SAM stdin preprocessed/
  """
}

process music_sort_forward {
  echo params.verbose
  publishDir "results/readthrough/sorted", mode: copy
  cpu 1
  input:
    file bams from preprocessed_forward_bams.collect()
  output:
    file "sorted/*" into sorted_forward_bams
  script:
  """
  mkdir to_sort
  mv ${bams} to_sort/
  mkdir sorted
  MUSIC -sort_reads ./to_sort/ sorted/
  """
}

process music_sort_reverse {
  echo params.verbose
  publishDir "results/readthrough/sorted", mode: copy
  cpu 1
  input:
    file bams from preprocessed_reverse_bams.collect()
  output:
    file "sorted/*" into sorted_reverse_bams
  script:
  """
  mkdir to_sort
  mv ${bams} to_sort/
  mkdir sorted
  MUSIC -sort_reads ./to_sort/ sorted/
  """
}

process music_deduplicat_forward {
  echo params.verbose
  publishDir "results/readthrough/deduplicated", mode: copy
  cpu 1
  input:
    file bams from sorted_forward_bams.collect()
  output:
    file "deduplicated/*" into deduplicated_forward_bams
  script:
  """
  mkdir to_deduplicate
  mv ${bams} to_deduplicate/
  mkdir deduplicated
  MUSIC -remove_duplicates ./to_deduplicate/ deduplicatede/
  """
}

process music_deduplicat_reverse {
  echo params.verbose
  publishDir "results/readthrough/deduplicated", mode: copy
  cpu 1
  input:
    file bams from sorted_reverse_bams.collect()
  output:
    file "deduplicated/*" into deduplicated_reverse_bams
  script:
  """
  mkdir to_deduplicate
  mv ${bams} to_deduplicate/
  mkdir deduplicated
  MUSIC -remove_duplicates ./to_deduplicate/ deduplicatede/
  """
}

deduplicated_reverse_bams.choice(
  wt_deduplicated_bams,
  mutants_deduplicated_bams){ a -> a =~ /.*_wt_.*/ ? 0 : 1 }

process music_computation {
  echo params.verbose
  publishDir "results/readthrough/deduplicated", mode: copy
  cpu 1
  input:
    file wt_bams from wt_deduplicated_bams.collect()
    file mutants_bams from mutants_deduplicated_bams.collect()
    file mapp from mappability.collect()
  output:
    file "peak_calling/*" into deduplicated_reverse_bams
  script:
  """
  mkdir mappability
  mv ${mapp} mappability/
  mkdir wt
  mv ${wt_bams} wt/
  mkdir mutants
  mv ${mutants_bams} mutants/
  MUSIC -get_per_win_p_vals_vs_FC -chip mutants/ -control wt/ \
    -l_win_step 50 -l_win_min 100 -l_win_max 1000
  MUSIC -get_multiscale_punctate_ERs \
    -chip mutants/ -control wt/ -mapp mappability/ \
    -begin_l 100 -end_l 10000 -step 1.2 \
    -l_mapp ${params.reads_size} -l_frag ${params.frag_size} -q_val 1 -l_p 0
  """
}
