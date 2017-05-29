#!/usr/bin/env nextflow
/*
* Copyright laurent modolo for the LBMC UMR 5239 ©.
* contributor(s) : laurent modolo (2017)
*
* laurent.modolo@ens-lyon.fr
*
* This software is a computer program whose purpose is to run the
* mapping and quantification steps of an NGS analysis of this project
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

rootDir = (baseDir =~ /(.*)src\/pipe/)[0][1]
results_path = rootDir + 'results'
quality_control_path = results_path + '/mapping'
fastqc_res_path = quality_control_path + '/bams'
multiqc_res_path = quality_control_path + '/counts'
src_path = rootDir + '/src'
params.name = "mapping and quantification analysis"
params.salmon = "/usr/bin/salmon"
if( !file(params.salmon).exists() ) exit 1, "salmon binary not found at: ${params.salmon}"
params.kallisto = "/usr/local/bin/kallisto"
if( !file(params.kallisto).exists() ) exit 1, "kallisto binary not found at: ${params.kallisto}"
params.paired = true

if(params.paired != true && params.paired != false){
   exit 1, "Invalid paired option: ${params.paired}. Valid options: 'true' or 'false'"
}

log.info params.name
log.info "============================================"
log.info "fastq files : ${params.fastq_files}"
log.info "paired files : ${params.paired}"
if (params.paired) {
  log.info "file names are expected to end in the format *_{1,2}.fastq*."
  log.info "or *_R{1,2}.fastq*."
  log.info "otherwise the pairs will not be paired for the analysis"
}
log.info "reference files : ${params.reference}"
log.info "salmon : ${params.salmon}"
log.info "kallisto : ${params.kallisto}"
log.info "results folder : ${results_path}"
log.info "\n"

if (params.paired){
    fastq_names = Channel.fromFilePairs( params.fastq_files, size:2)
    .ifEmpty { exit 1, "Cannot find any fastq files matching: ${params.fastq_files}" }
} else {
    fastq_names = Channel.fromPath( params.fastq_files )
    .ifEmpty { exit 1, "Cannot find any fastq files matching: ${params.fastq_files}" }
}

process get_file_name_fastq {
  tag "${tagname}"
  input:
    val file_name from fastq_names
  output:
    file "*${file_name_root}" into dated_fastq_names
  when:
    if (params.paired) {
      (file_name[1][0] =~ /^.*\.fastq$/ || file_name[1][0] =~ /^.*\.fastq\.gz$/) && (file_name[1][1] =~ /^.*\.fastq$/ || file_name[1][1] =~ /^.*\.fastq\.gz$/)
    } else {
      file_name =~ /^.*\.fastq$/ || file_name =~ /^.*\.fastq\.gz$/
    }
  script:
  if (params.paired) {
    tagname = file_name[0]
    file_name_root = file_name[0] + "*"
    """
    ${src_path}/func/file_handle.py -f ${file_name[1][0]} ${file_name[1][1]} -c -e | \
    awk '{system("ln -s "\$0" ."); print(\$0)}'
    """
  } else {
    tagname = file(file_name).baseName
    file_name_root = file_name.name
    """
    ${src_path}/func/file_handle.py -f ${file_name} -c -e | \
    awk '{system("ln -s "\$0" ."); print(\$0)}'
    """
  }
}

reference_names = Channel.fromPath( params.reference_files)
.ifEmpty { exit 1, "Cannot find any reference file matching: ${params.reference_files}" }
process get_file_name_fasta {
  tag "${tagname}"
  input:
    val file_name from reference_names
  output:
    file "*${file_name_root}" into dated_reference_names
  when:
    file_name =~ /^.*\.fasta$/ || file_name =~ /^.*\.fasta\.index$/
  script:
    tagname = file(file_name).baseName
    file_name_root = file_name.name
    """
    ${src_path}/func/file_handle.py -f ${file_name} -c -e | \
    awk '{system("ln -s "\$0" ."); print(\$0)}'
    """
}
