#!/usr/bin/env nextflow
/*
* Copyright laurent modolo for the LBMC UMR 5239 ©.
* contributor(s) : laurent modolo (2017)
*
* laurent.modolo@ens-lyon.fr
*
* This software is a computer program whose purpose is to test the
* src/func/file_handle.py program of this project
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

results_path = baseDir + '/../../results'
quality_control_path = results_path + '/quality_control'
fastqc_res_path = quality_control_path + '/fastqc'
multiqc_res_path = quality_control_path + '/multiqc'
adaptor_removal_res_path = quality_control_path + '/adaptor_removed'
trimming_res_path = quality_control_path + '/trimming'
src_path = baseDir + '/../../src'
params.name = "quality control analysis"
params.fastqc = "/usr/bin/fastqc"
if( !file(params.fastqc).exists() ) exit 1, "fastqc binary not found at: ${params.fastqc}"
params.multiqc = "/usr/bin/multiqc"
if( !file(params.multiqc).exists() ) exit 1, "multiqc binary not found at: ${params.multiqc}"
paired = false
params.adaptor_removal = "cutadapt"
params.adaptor_sequence = "-a AGATCGGAAGAG -g CTCTTCCGATCT"
if(paired){
  params.adaptor_sequence = "-a AGATCGGAAGAG -g CTCTTCCGATCT -A AGATCGGAAGAG -G CTCTTCCGATCT"
}
params.trimmer = "UrQt"
params.quality_threshold = 20
params.urqt = "/usr/bin/UrQt"
params.trimmomatic = "/usr/bin/trimmomatic"
params.cutadapt = "/usr/bin/cutadapt"

log.info params.name
log.info "============================================"
log.info "fastqc : ${params.fastqc}"
log.info "multiqc : ${params.multiqc}"
log.info "query : ${params.fastq_files}"
log.info "results folder : ${results_path}"
log.info "fastq_files : ${params.fastq_files}"
log.info "trimmer : ${params.trimmer}"
if (params.trimmer == 'cutadapt') {
  log.info "cutadapt path : ${params.cutadapt}"
}else{
  log.info "UrQt path : ${params.urqt}"
}
log.info "\n"

Channel
  .from( params.fastq_files )
  .set{ file_names }
Channel
  .create()
  .set{ dated_file_names }
Channel
  .create()
  .set{ dated_fastqc_files }
Channel
  .create()
  .set{ multiqc_report }

process get_file_name {
  input:
    val file_name from file_names
  output:
    stdout dated_file_name into dated_file_names
  script:
  """
  echo -e \$(${src_path}/func/file_handle.py -f $baseDir/../../${file_name} -c -e)
  """
}

dated_file_names.splitCsv().map{ n -> tuple(n, n) }.into{ fastqc_input; adaptor_rm_input }

process fastqc {
  tag "${file(name[0]).name}"
  publishDir "${fastqc_res_path}", mode: 'copy'
  input:
    set val(name), val(file_name) from fastqc_input
  output:
    file "*_fastqc.{zip,html}" into fastqc_output
  when:
    file(file_name[0]).name =~ /^.*\.fastq$/ || file(file_name[0]).name =~ /^.*\.fastq\.gz$/
  script:
  """
    ${params.fastqc} --quiet --outdir ./ ${file_name[0]}
    ${src_path}/func/file_handle.py -f *.html -r
    ${src_path}/func/file_handle.py -f *.zip -r
  """
}

process adaptor_removal {
  tag "${file(name[0]).name}"
  publishDir "${adaptor_removal_res_path}", mode: 'copy'
  input:
    set val(name), val(file_name) from adaptor_rm_input
  output:
    file "*.fastq.gz" into adaptor_rm_output
  when:
    file(file_name[0]).name =~ /^.*\.fastq$/ || file(file_name[0]).name =~ /^.*\.fastq\.gz$/
  script:
  if (params.adaptor_removal == "cutadapt") {
  """
    ${params.cutadapt} ${params.adaptor_sequence} ${file_name[0]} -o ${file(file_name[0]).baseName}.fastq.gz
    ${src_path}/func/file_handle.py -f *.fastq.gz -r
  """
  }
}

adaptor_rm_output.map{ n -> tuple(n, n) }.into{ trimming_input }

process trimming {
  tag "${name.name}"
  publishDir "${trimming_res_path}", mode: 'copy'
  input:
    set file(name), file(file_name) from trimming_input
  output:
    file "*.trimmed.fastq.gz" into trimming_output
    file "*_report.txt" into trimming_log
  when:
    file_name.name =~ /^.*\.fastq$/ || file_name.name =~ /^.*\.fastq\.gz$/
  script:
  if (params.trimmer == "cutadapt") {
  """
    ${params.cutadapt} -q ${params.quality_threshold},${params.quality_threshold} ${file_name} -o ${file_name.baseName}.trimmed.fastq.gz > ${file_name.baseName}_report.txt
    ${src_path}/func/file_handle.py -f *.trimmed.fastq.gz -r
  """
  }else{
  """
    ${params.urqt} --t ${params.quality_threshold} --gz --in ${file_name} --out ${file_name.baseName}.trimmed.fastq.gz > ${file_name.baseName}_report.txt
    ${src_path}/func/file_handle.py -f *.trimmed.fastq.gz -r
  """
  }
}

trimming_output.map{ n -> tuple(n, n) }.into{ fastqc_trimmed_input }

process fastqc_trimmed {
  tag "${name.name}"
  publishDir "${fastqc_res_path}", mode: 'copy'
  input:
    set file(name), file(file_name) from fastqc_trimmed_input
  output:
    file "*_fastqc.{zip,html}" into fastqc_trimmed_output
  when:
    file_name.name =~ /^.*\.fastq$/ || file_name.name =~ /^.*\.fastq\.gz$/
  script:
  """
    ${params.fastqc} --quiet --outdir ./ ${file_name}
    ${src_path}/func/file_handle.py -f *.html -r
    ${src_path}/func/file_handle.py -f *.zip -r
  """
}

process multiqc {
    publishDir "${multiqc_res_path}", mode: 'copy'
    input:
      file fastqc_results from fastqc_output.collect()
      file fastqc_trimmed_results from fastqc_trimmed_output.collect()
    output:
      file "*multiqc_{report,data}*" into multiqc_report
    script:
    """
    ${params.multiqc} -f .
    ${src_path}/func/file_handle.py -f multiqc_report.html -r
    ${src_path}/func/file_handle.py -f multiqc_data -r
    """
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
