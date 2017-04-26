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

results_path = baseDir + '/../../results/'
src_path = baseDir + '/../../src/'
params.name = "quality control analysis"
params.fastqc = "/usr/bin/fastqc"
fastqc_res_path = results_path+'/quality_control/fastqc/'
if( !file(params.fastqc).exists() ) exit 1, "fastqc binary not found at: ${params.fastqc}"
params.multiqc = "/usr/bin/multiqc"
multiqc_res_path = results_path+'/quality_control/multiqc/'
if( !file(params.multiqc).exists() ) exit 1, "multiqc binary not found at: ${params.multiqc}"

log.info params.name
log.info "============================================"
log.info "fastqc : ${params.fastqc}"
log.info "multiqc : ${params.multiqc}"
log.info "query : ${params.fastq_files}"
log.info "results folder : ${results_path}"
fastq_files = params.fastq_files.tokenize(' ')
log.info "fastq_files : ${fastq_files}"
log.info "\n"

/*
We get the dated names of the files and we send them to:
fastqc
trimming programme
*/
dated_fastq_files_fastqc = Channel.create()
dated_fastq_files_trimming = Channel.create()
dated_fastq_files_names = Channel.create()
dated_fastq_files_names.flatMap{ n -> n.split("\n") }
                       .map{ n -> file(n) }
                       .into( dated_fastq_files_fastqc )
process get_file_name {
  input:
  val fastq_file from fastq_files
  output:
  stdout dated_fastq_file into dated_fastq_files_names
  script:
  """
  echo -e \$(${src_path}/func/file_handle.py --file $baseDir/../../${fastq_file} --check --escape)
  """
}

/*
We run fastqc on a list for fastq files
*/
fastqc_results = Channel.create()
process fastqc {
  publishDir fastqc_res_path, mode: 'copy'
  input:
  file dated_fastq_file from dated_fastq_files_fastqc
  output:
  file "*_fastqc.{zip,html}" into fastqc_results
  when:
  dated_fastq_file.extension =~ /^fastq/ || dated_fastq_file_to_qc.extension =~ /^fastq\.gz/
  script:
  """
    ${params.fastqc} --quiet $dated_fastq_file
    ${src_path}/func/file_handle.py --file ${dated_fastq_file.baseName}_fastqc.* --redate
  """
}

/*
We run multiqc on the results of fastqc
*/
process multiqc {
    publishDir "${multiqc_res_path}", mode: 'copy'
    input:
    file ("$fastqc_res_path/*") from fastqc_results.flatten().toList()
    output:
    file '*multiqc_report.html' into multiqc_html
    file '*multiqc_data' into multiqc_data
    script:
    """
    python2 ${params.multiqc} -f .
    """
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
