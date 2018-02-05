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
params.fastq = ""
params.mean_size = 200
params.sd_size = 200

fastq_files = Channel
  .fromPath(params.fastq)
  .ifEmpty {
    exit 1, "Cannot find any fastq file matching: ${params.fastq}"
  }

fasta_files = Channel
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
    file annotation from annotation_forward.collect()
  output:
    file "*.gff" into gff3_forward
  script:
    """
    cat ${annotation} > annotation_merge.bed

    grep "transcript" annotation_merge.bed | \
      grep -v "transcript_RT" | \
      awk '{if(\$2 < \$3) {print \$0}}' > \
      annotation_merge_T.bed

    grep "transcript_RT" annotation_merge | \
      awk '{if(\$2 < \$3) {print \$0}}' > \
      annotation_merge_RT.bed

    cat annotation_merge_T.bed annotation_merge_RT.bed | \
      bedtools sort -i stdin | \
      gt-bed_to_gff3 -o forward.gff -i stdin

    file_handle.py -f *.gff
    """
}

process merge_annotation_reverse {
  echo params.verbose
  publishDir "results/readthrough/transcript/", mode: 'copy'
  cpus 1
  input:
    file annotation from annotation_reverse.collect()
  output:
    file "*.gff" into gff3_reverse
  script:
    """
    cat ${annotation} > annotation_merge.bed

    grep "transcript" annotation_merge.bed | \
      grep -v "transcript_RT" | \
      awk '{if(\$2 < \$3) {print \$0}}' > \
      annotation_merge_T.bed

    grep "transcript_RT" annotation_merge | \
      awk '{if(\$2 < \$3) {print \$0}}' > \
      annotation_merge_RT.bed

    cat annotation_merge_T.bed annotation_merge_RT.bed | \
      bedtools sort -i stdin | \
      gt-bed_to_gff3 -o reverse.gff -i stdin

    file_handle.py -f *.gff
    """
}
