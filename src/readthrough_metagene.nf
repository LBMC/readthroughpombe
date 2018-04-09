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
params.annotation = ""
params.bam = ""

bam_files = Channel
  .fromPath(params.bam)
  .ifEmpty {
    exit 1, "Cannot find any bam file matching: ${params.bam}"
  }

annotation_files = Channel
  .fromPath(params.annotation)
  .ifEmpty {
    exit 1, "Cannot find any annotation file matching: ${params.annotation}"
  }

process merge_annotation {
  echo params.verbose
  publishDir "results/readthrough/metagene/bed/", mode: 'copy'
  cpus 1
  input:
    file annotation from annotation_files.collect()
  output:
    file "*RT_annotation.bed" into rt_bed
  script:
    """
    cat ${annotation} | \
      awk '{if(\$2 < \$3) {print \$0}}' | \
      bedtools sort -i stdin > RT_annotation.bed

    file_handle.py -f RT_annotation.bed
    """
}

bam_files.groupBy {
  str -> (str =~ /^.*\d{5}_\d{2}_\d{2}_[a-zA-Z0-9]+_(.*)_.*$/)[0][1]
}
.flatMap().
map{
  it -> it.getValue()
}
.set{bams_grouped}

process merge_bams {
  echo params.verbose
  cpus 1
  input:
    file bams from bams_grouped
  output:
    file "*.bam*" into bams_merged
  script:
    condition = (bams[0] =~ /^\d{4}_\d{2}_\d{2}_(.*)_[a-zA-Z0-9]+\..*/)[0][1]
    """
    samtools merge ${condition}.bam ${bams}
    file_handle.py -f ${condition}.bam
    """
}

process sort_bams {
  echo params.verbose
  publishDir "results/readthrough/metagene/bams/", mode: 'copy', pattern: "_sort.bam"
  cpus 4
  input:
    file bams from bams_merged
  output:
    file "*.bam*" into bams_sorted
  script:
    file_name = (bams =~ /^(.*)\.bam/)[0][1]
    """
    samtools sort -@ ${task.cpus} -o ${file_name}_sort.bam ${bams}
    file_handle.py -f ${file_name}_sort.bam
    """
}

process index_bams {
  echo params.verbose
  publishDir "results/readthrough/metagene/bams/", mode: 'copy', pattern: ".bam.bai"
  cpus 1
  input:
    file bams from bams_sorted
  output:
    file "*.bam*" into ibam_files
  script:
    """
    samtools index ${bams}
    mv ${bams} ${bams}d
    file_handle.py -f ${bams}*
    """
}

process compute_bigwig {
  echo params.verbose
  publishDir "results/readthrough/metagene/bigwig/", mode: 'copy'
  cpus 11
  input:
    set file(bais), file(bams) from ibam_files
  output:
    file "*.bw" into rt_bw
  script:
    bams = (bams =~ /^(.*)d$/)[0][1]
    """
    mv ${bams}d ${bams}
    bamCoverage --bam ${bams} --outFileFormat bigwig -o ${bams}.bw\
      --binSize 10 -p ${task.cpus} --normalizeUsing BPM
    file_handle.py -f ${bams}.bw
    """
}

process compute_matrix {
  echo params.verbose
  publishDir "results/readthrough/metagene/matrix/", mode: 'copy'
  cpus 11
  input:
    file bw from rt_bw.collect()
    file bed from rt_bed
  output:
    file "*.mat" into rt_mat
  script:
    """
    computeMatrix reference-point -S ${bw} -R ${bed} --referencePoint TSS -b 200 -a 200 -o metagene.mat -p ${task.cpus}
    file_handle.py -f metagene.mat
    """
}

process compute_plot {
  echo params.verbose
  publishDir "results/readthrough/metagene/img/", mode: 'copy'
  cpus 1
  input:
    file mat from rt_mat
  output:
    file "*.pdf" into rt_pdf
  script:
    """
    plotProfile -m ${mat} --perGroup -o ${mat}.pdf --plotTitle "metagene"
    plotHeatmap -m ${mat} --kmean 2 -out ${mat}_hm.pdf
    file_handle.py -f ${mat}.pdf ${mat}_hm.pdf
    """
}

