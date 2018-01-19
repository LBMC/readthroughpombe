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
params.fastq = ""

fastq_files = Channel.fromPath(params.fastq)
  .ifEmpty { exit 1, "Cannot find any fastq files matching: ${params.fastq}" }

process reverse_fastq {
  echo params.verbose
  publishDir "results/reversecomplement/", mode: 'copy'
  cpus = 1
  input:
    file fastq from fastq_files
  output:
    file "*_rev.fastq.gz" into reverse_fastq
  script:
    """
    seqkit seq -r -p ${fastq} | gzip -c  > ${fastq}_rev
    find . -name "*_rev" | sed 's/\\.fastq\\.gz_rev//g' | \
    awk '{system("mv "\$0".fastq.gz_rev "\$0"_rev.fastq.gz")}'
    file_handle.py -f *.fastq.gz
    """
}
