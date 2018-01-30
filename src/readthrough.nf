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
params.name = ""
params.bams = ""
// params.bams_mutant = ""
// params.bams_wt = ""
params.annotation = ""
params.index = ""
params.genome = ""
params.reads_size = 0
params.frag_size = 200
params.begin_l = 50
params.end_l = 500

bam_files = Channel
    .fromPath(params.bams, glob: true)
    .splitCsv( header: true, strip: true )
    .flatMap{
       it -> [
        [file(it['mutant']), it['mutant'], "mutant"],
        [file(it['wt']), it['mutant'], 'wt']
      ]
    }

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
  publishDir "results/readthrough/bams/" + params.name, mode: 'copy'
  cpus 4
  input:
    set file(bams), val(pair), val(type) from bam_files
  output:
    set file("*_sorted.bam"), val(pair), val(type) into sorted_bams
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
  publishDir "results/readthrough/bams/" + params.name, mode: 'copy'
  cpus 2
  input:
    set file(bams), val(pair), val(type) from sorted_bams
  output:
    set file("forward/*_forward.bam*"), val(pair), val(type) into forward_bams
    set file("reverse/*_reverse.bam*"), val(pair), val(type) into reverse_bams
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

process gff_to_bed {
  echo params.verbose
  publishDir "results/readthrough/bams/" + params.name, mode: 'copy'
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
  file_handle.py -f *_reverse.bed *_forward.bed
  """
}

genome_file.into{
  genome_file_mappability;
  genome_file_to_annotate;
  genome_file_transcript_forward;
  genome_file_transcript_reverse
}

annotation_forward.into{
  annotation_forward_to_merge;
  annotation_forward_to_negate
}

annotation_reverse.into{
  annotation_reverse_to_merge;
  annotation_reverse_to_negate
}

process build_genome {
  echo params.verbose
  publishDir "results/readthrough/bams/" + params.name + "/", mode: 'copy'
  cpus 4
  input:
    file genome from genome_file_to_annotate
  output:
    file "*.genome" into genome_file_annotation
  script:
  """
  samtools faidx ${genome}
  awk -v OFS='\t' '{print \$1,\$2}' ${genome}.fai | \
  sort -k1 > ${genome}.genome
  file_handle.py -f *.genome
  """
}

genome_file_annotation.into{
  genome_file_annotation_forward;
  genome_file_annotation_reverse;
  genome_file_peak_merge_forward;
  genome_file_peak_merge_reverse
}

process negative_forward {
  echo params.verbose
  publishDir "results/readthrough/bams/" + params.name + "/forward", mode: 'copy'
  cpus 4
  input:
    file genome from genome_file_annotation_forward
    file annotation from annotation_forward_to_negate
  output:
    file "*_negative.bed" into annotation_negative_forward
  script:
  """
  bedtools complement -i ${annotation} -g ${genome} > ${annotation}_negative.bed
  find . -name "*_negative.bed" | sed 's/\\.bed_negative\\.bed//g' | \
  awk '{system("mv "\$0".bed_negative.bed "\$0"_negative.bed")}'
  file_handle.py -f *_negative.bed
  """
}

process negative_reverse {
  echo params.verbose
  publishDir "results/readthrough/bams/" + params.name + "/reverse", mode: 'copy'
  cpus 4
  input:
    file genome from genome_file_annotation_reverse
    file annotation from annotation_reverse_to_negate
  output:
    file "*_negative.bed" into annotation_negative_reverse
  script:
  """
  bedtools complement -i ${annotation} -g ${genome} > ${annotation}_negative.bed
  find . -name "*_negative.bed" | sed 's/\\.bed_negative\\.bed//g' | \
  awk '{system("mv "\$0".bed_negative.bed "\$0"_negative.bed")}'
  file_handle.py -f *_negative.bed
  """
}

process filter_forward {
  echo params.verbose
  publishDir "results/readthrough/bams/" + params.name + "/forward", mode: 'copy'
  cpus 4
  input:
    set file(bams), val(pair), val(type) from forward_bams
    file annotation from annotation_negative_forward.first()
  output:
    set file("*_noannot.bam*"), val(pair), val(type) into filtered_forward_bams
  script:
  """
    samtools view -@ ${task.cpus} -hb ${bams} -L ${annotation} > ${bams}_noannot.bam
    samtools index ${bams}_noannot.bam
    find . -name "*_noannot.bam" | sed 's/\\.bam_noannot\\.bam//g' | \
    awk '{system("mv "\$0".bam_noannot.bam "\$0"_noannot.bam")}'
    find . -name "*_noannot.bam.bai" | sed 's/\\.bam_noannot\\.bam.bai//g' | \
    awk '{system("mv "\$0".bam_noannot.bam.bai "\$0"_noannot.bam.bai")}'
    file_handle.py -f *_noannot.bam*
  """
}

process filter_reverse {
  echo params.verbose
  publishDir "results/readthrough/bams/" + params.name + "/reverse", mode: 'copy'
  cpus 4
  input:
    set file(bams), val(pair), val(type) from reverse_bams
    file annotation from annotation_negative_reverse.first()
  output:
    set file("*_noannot.bam*"), val(pair), val(type) into filtered_reverse_bams
  script:
  """
    samtools view -@ ${task.cpus} -hb ${bams} -L ${annotation} > ${bams}_noannot.bam
    samtools index ${bams}_noannot.bam
    find . -name "*_noannot.bam" | sed 's/\\.bam_noannot\\.bam//g' | \
    awk '{system("mv "\$0".bam_noannot.bam "\$0"_noannot.bam")}'
    find . -name "*_noannot.bam.bai" | sed 's/\\.bam_noannot\\.bam.bai//g' | \
    awk '{system("mv "\$0".bam_noannot.bam.bai "\$0"_noannot.bam.bai")}'
    file_handle.py -f *_noannot.bam*
  """
}

process compute_mappability {
  echo params.verbose
  cpus 4
  input:
    file index from index_file.collect()
    file genome from genome_file_mappability
  output:
    file "*.bin" into mappability
  script:
  """
  generate_multimappability_signal.csh ${genome} ${params.reads_size} ./
  bash temp_map_reads.csh
  bash temp_process_mapping.csh
  """
}

filtered_forward_bams.into{
  filtered_forward_bams_preprocess;
  filtered_forward_bams_sort;
  filtered_forward_bams_deduplicate;
  filtered_forward_bams_peak_calling
}

filtered_reverse_bams.into{
  filtered_reverse_bams_preprocess;
  filtered_reverse_bams_sort;
  filtered_reverse_bams_deduplicate;
  filtered_reverse_bams_peak_calling
}

process music_preprocess_forward {
  echo params.verbose
  cpus 1
  input:
    set file(bams), val(pair), val(type) from filtered_forward_bams_preprocess
  output:
    file "preprocessed/*" into preprocessed_forward
  script:
  """
  mkdir preprocessed
  samtools view *.bam | \
  MUSIC -preprocess SAM stdin preprocessed/
  find preprocessed/ | awk '{system("mv "\$0" "\$0"_tmp")}'
  """
}

process music_preprocess_reverse {
  echo params.verbose
  cpus 1
  input:
    set file(bams), val(pair), val(type) from filtered_reverse_bams_preprocess
  output:
    file "preprocessed/*" into preprocessed_reverse
  script:
  """
  mkdir preprocessed
  samtools view *.bam | \
  MUSIC -preprocess SAM stdin preprocessed/
  find preprocessed/ | awk '{system("mv "\$0" "\$0"_tmp")}'
  """
}

process music_sort_forward {
  echo params.verbose
  cpus 1
  input:
    set file(bams), val(pair), val(type) from filtered_forward_bams_sort
    file preprocess from preprocessed_forward
  output:
    file "sorted/*" into sorted_forward
  script:
  """
  mkdir to_sort
  mv ${preprocess} to_sort/
  find to_sort | sed 's/_tmp//g' |awk '{system("mv "\$0"_tmp "\$0)}'
  mv ${bams} to_sort/
  mkdir sorted
  MUSIC -sort_reads to_sort/ sorted/
  find sorted/ | awk '{system("mv "\$0" "\$0"_tmp2")}'
  """
}

process music_sort_reverse {
  echo params.verbose
  cpus 1
  input:
    set file(bams), val(pair), val(type) from filtered_reverse_bams_sort
    file preprocess from preprocessed_reverse
  output:
    file "sorted/*" into sorted_reverse
  script:
  """
  mkdir to_sort
  mv ${preprocess} to_sort/
  find to_sort | sed 's/_tmp//g' |awk '{system("mv "\$0"_tmp "\$0)}'
  mv ${bams} to_sort/
  mkdir sorted
  MUSIC -sort_reads ./to_sort/ sorted/
  find sorted/ | awk '{system("mv "\$0" "\$0"_tmp2")}'
  """
}

process music_deduplicat_forward {
  echo params.verbose
  cpus 1
  input:
    set file(bams), val(pair), val(type) from filtered_forward_bams_deduplicate
    file sorted from sorted_forward
  output:
    file "deduplicated/*" into deduplicated_forward
  script:
  """
  mkdir to_deduplicate
  mv ${sorted} to_deduplicate/
  find to_deduplicate/* | sed 's/_tmp2//g' |awk '{system("mv "\$0"_tmp2 "\$0)}'
  mv ${bams} to_deduplicate/
  mkdir deduplicated
  MUSIC -remove_duplicates ./to_deduplicate/ 2 deduplicated/
  if [ -f to_deduplicate/*_wt_*.bam ]; then
    find deduplicated/* | awk '{system("mv "\$0" "\$0"_wt")}'
  fi
  if [ ! -f to_deduplicate/*_wt_*.bam ]; then
    find deduplicated/* | awk '{system("mv "\$0" "\$0"_mutant")}'
  fi
  """
}

process music_deduplicat_reverse {
  echo params.verbose
  cpus 1
  input:
    set file(bams), val(pair), val(type) from filtered_reverse_bams_deduplicate
    file sorted from sorted_reverse
  output:
    file "deduplicated/*" into deduplicated_reverse
  script:
  """
  mkdir to_deduplicate
  mv ${sorted} to_deduplicate/
  find to_deduplicate/* | sed 's/_tmp2//g' |awk '{system("mv "\$0"_tmp2 "\$0)}'
  mv ${bams} to_deduplicate/
  mkdir deduplicated
  MUSIC -remove_duplicates ./to_deduplicate/ 2 deduplicated/
  if [ -f to_deduplicate/*_wt_*.bam ]; then
    find deduplicated/* | awk '{system("mv "\$0" "\$0"_wt")}'
  fi
  if [ ! -f to_deduplicate/*_wt_*.bam ]; then
    find deduplicated/* | awk '{system("mv "\$0" "\$0"_mutant")}'
  fi
  """
}

mappability.collect().into{
  mappability_forward;
  mappability_reverse
}


filtered_forward_bams_peak_calling
  .groupTuple(by: 1)
  .map{
    it -> [it[0][0], it[0][1]]
  }.set{mutant_wt_forward_bams}

filtered_reverse_bams_peak_calling
  .groupTuple(by: 1)
  .map{
    it -> [it[0][0], it[0][1]]
  }.set{mutant_wt_reverse_bams}

wt_deduplicated_forward =  Channel.create()
mutant_deduplicated_forward =  Channel.create()

deduplicated_forward.choice(
  wt_deduplicated_forward,
  mutant_deduplicated_forward){ a -> a =~ /.*_wt/ ? 0 : 1 }

wt_deduplicated_reverse =  Channel.create()
mutant_deduplicated_reverse =  Channel.create()

deduplicated_reverse.choice(
  wt_deduplicated_reverse,
  mutant_deduplicated_reverse){ a -> a =~ /.*_wt/ ? 0 : 1 }

process music_forward_computation {
  echo params.verbose
  publishDir "results/readthrough/peak_calling/" + params.name + "/forward", mode: 'copy'
  memory '30GB'
  cpus 1
  input:
    set file(mutant_bams), file(wt_bams) from mutant_wt_forward_bams
    file wt_deduplicated from wt_deduplicated_forward
    file mutant_deduplicated from mutant_deduplicated_forward
    file mapp from mappability_forward.first()
  output:
    file "*" into music_output_forward
    file "*_ERs*.bed" into peaks_forward
  script:
  """
  mkdir mappability
  mv ${mapp} mappability/

  mkdir wt
  mv ${wt_deduplicated} wt/
  find wt/ | sed 's/_wt//g' |awk '{system("mv "\$0"_wt "\$0)}'
  mv ${wt_bams} wt/

  mkdir mutant
  mv ${mutant_deduplicated} mutant/
  find mutant/ | sed 's/_mutant//g' |awk '{system("mv "\$0"_mutant "\$0)}'
  mv ${mutant_bams} mutant/

  MUSIC -get_per_win_p_vals_vs_FC -chip mutant/ -control wt/ \
    -l_win_step 50 -l_win_min 50 -l_win_max 5000
  MUSIC -get_multiscale_punctate_ERs \
    -chip mutant/ -control wt/ -mapp mappability/ \
    -begin_l ${params.begin_l} -end_l ${params.end_l} -step 1.1 \
    -l_mapp ${params.reads_size} -l_frag ${params.frag_size} -q_val 1 -l_p 0

  find mutant/*.bam |\
    sed 's/\\.bam//g' |\
    sed 's/mutant\\///g' |\
    sed 's/\\(.*\\)_rev*/\\1/g' |\
    sed 's/\\(.*\\)_trim*/\\1/g' |\
    awk '{system("mkdir "\$0"; rm -Rf wt mutant; mv * "\$0"; cp "\$0"/ERs*.bed ./")}'
  file_handle.py -f *
  """
}

process music_reverse_computation {
  echo params.verbose
  publishDir "results/readthrough/peak_calling/" + params.name + "/reverse", mode: 'copy'
  cpus 1
  memory '30GB'
  input:
    set file(mutant_bams), file(wt_bams) from mutant_wt_reverse_bams
    file wt_deduplicated from wt_deduplicated_reverse
    file mutant_deduplicated from mutant_deduplicated_reverse
    file mapp from mappability_reverse.first()
  output:
    file "*" into music_output_reverse
    file "*_ERs*.bed" into peaks_reverse
  script:
  """
  mkdir mappability
  mv ${mapp} mappability/

  mkdir wt
  mv ${wt_deduplicated} wt/
  find wt/ | sed 's/_wt//g' |awk '{system("mv "\$0"_wt "\$0)}'
  mv ${wt_bams} wt/

  mkdir mutant
  mv ${mutant_deduplicated} mutant/
  find mutant/ | sed 's/_mutant//g' |awk '{system("mv "\$0"_mutant "\$0)}'
  mv ${mutant_bams} mutant/

  MUSIC -get_per_win_p_vals_vs_FC -chip mutant/ -control wt/ \
    -l_win_step 50 -l_win_min 50 -l_win_max 5000
  MUSIC -get_multiscale_punctate_ERs \
    -chip mutant/ -control wt/ -mapp mappability/ \
    -begin_l ${params.begin_l} -end_l ${params.end_l} -step 1.1 \
    -l_mapp ${params.reads_size} -l_frag ${params.frag_size} -q_val 1 -l_p 0

  find mutant/*.bam |\
    sed 's/\\.bam//g' |\
    sed 's/mutant\\///g' |\
    sed 's/\\(.*\\)_rev*/\\1/g' |\
    sed 's/\\(.*\\)_trim*/\\1/g' |\
    awk '{system("mkdir "\$0"; rm -Rf wt mutant; mv * "\$0"; cp "\$0"/ERs*.bed ./")}'
  file_handle.py -f *
  """
}

process bed_merge_forward {
  echo params.verbose
  publishDir "results/readthrough/peak_calling/", mode: 'copy'
  echo true
  input:
    file bed from peaks_forward.collect()
    file genome from genome_file_peak_merge_forward
  output:
    file "*_forward_full.bed" into peaks_final_forward_full
    file "*_forward.bed" into peaks_final_forward
  script:
  """
  cat ${bed} > concatenate_forward_m.bed
  bedtools sort -i concatenate_forward_m.bed | \
  awk -v OFS='\t' '{print \$1, \$2, \$3, "-"}' > concatenate_forward_s.bed

  bedtools genomecov -bg -i concatenate_forward_s.bed -g ${genome} | \
  awk -v OFS='\t' '{if(\$4 >= ${bed.size()/2}) {print \$0}}' \
  > concatenate_forward_cov.bed

  bedtools merge -d ${params.reads_size} -i concatenate_forward_cov.bed \
  > concatenate_forward_u.bed
  awk -v OFS='\t' '{ print \$0, ".", ".", "+", ".", "RT", ".", "."}' \
  concatenate_forward_u.bed > ${params.name}_peak_calling_forward.bed

  cp concatenate_forward_s.bed ${params.name}_peak_calling_forward_full_${params.begin_l}_${params.end_l}.bed
  file_handle.py -f *_forward.bed *_forward_full_${params.begin_l}_${params.end_l}.bed
  """
}

process bed_merge_reverse {
  echo params.verbose
  publishDir "results/readthrough/peak_calling/", mode: 'copy'
  echo true
  input:
    file bed from peaks_reverse.collect()
    file genome from genome_file_peak_merge_reverse
  output:
    file "*_reverse_full.bed" into peaks_final_reverse_full
    file "*_reverse.bed" into peaks_final_reverse
  script:
  """
  cat ${bed} > concatenate_reverse_m.bed
  bedtools sort -i concatenate_reverse_m.bed | \
  awk -v OFS='\t' '{print \$1, \$2, \$3, "-"}' > concatenate_reverse_s.bed

  bedtools genomecov -bg -i concatenate_reverse_s.bed -g ${genome} | \
  awk -v OFS='\t' '{if(\$4 >= ${bed.size()/2}) {print \$0}}' \
  > concatenate_reverse_cov.bed

  bedtools merge -d ${params.reads_size} -i concatenate_reverse_cov.bed \
  > concatenate_reverse_u.bed
  awk -v OFS='\t' '{ print \$0, ".", ".", "-", ".", "RT", ".", "."}' \
  concatenate_reverse_u.bed > ${params.name}_peak_calling_reverse.bed

  cp concatenate_reverse_s.bed ${params.name}_peak_calling_reverse_full_${params.begin_l}_${params.end_l}.bed
  file_handle.py -f *_reverse.bed *_reverse_full_${params.begin_l}_${params.end_l}.bed
  """
}

process peak_merge_forward {
  echo params.verbose
  publishDir "results/readthrough/peak_calling/", mode: 'copy'
  cpus 1
  input:
    file peaks from peaks_final_forward
    file annotation from annotation_forward_to_merge
  output:
    file "*RT_forward.bed" into rt_forward
  script:
  """
  bedtools substract -a ${peaks} -b ${annotation} > peaks_cleaned.bed
  cat peaks_cleaned.bed ${annotation} > to_merge_forward.bed
  bedtools sort -i to_merge_forward.bed > to_merge_forward_s.bed

  tac to_merge_forward_s.bed | \
  awk -v OFS='\t' '{
    if (\$8 == "RT") {
      rt_start = \$2;
      rt_stop = \$3;
    } else {
      if (\$8 == "transcript" && rt_start != "" && rt_start <= (\$3 + 100) ) {
        print \$0;
        \$3 = rt_stop;
        rt_start = "";
        rt_stop = "";
        \$8 = "transcript_RT";
        print \$0;
      } else {
        print \$0;
      }
    }
  }' > ${params.name}_RT_forward.bed

  file_handle.py -f ${params.name}_RT_forward.bed
  """
}

process peak_merge_reverse {
  echo params.verbose
  publishDir "results/readthrough/peak_calling/", mode: 'copy'
  cpus 1
  input:
    file peaks from peaks_final_reverse
    file annotation from annotation_reverse_to_merge
  output:
    file "*RT_reverse.bed" into rt_reverse
  script:
  """
  bedtools substract -a ${peaks} -b ${annotation} > peaks_cleaned.bed
  cat peaks_cleaned.bed ${annotation} > to_merge_reverse.bed
  bedtools sort -i to_merge_reverse.bed > to_merge_reverse_s.bed

  awk -v OFS='\t' '{
    if (\$8 == "RT") {
      rt_start = \$2;
      rt_stop = \$3;
    } else {
      if (\$8 == "transcript" && rt_start != "" && rt_stop <= (\$3 - 100) ) {
        print \$0
        \$2 = rt_start;
        rt_start = "";
        rt_stop = "";
        \$8 = "transcript_RT";
        print \$0;
      } else {
        print \$0;
      }
    }
  }' to_merge_reverse_s.bed > ${params.name}_RT_reverse.bed

  file_handle.py -f ${params.name}_RT_reverse.bed
  """
}
