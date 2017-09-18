#!/usr/bin/env nextflow
/*
* Copyright laurent modolo for the LBMC UMR 5239 ©.
* contributor(s) : laurent modolo (2017)
*
* laurent.modolo@ens-lyon.fr
*
* This software is a computer program whose purpose is to run the
* quality control steps of an NGS analysis of this project
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
quality_control_path = results_path + '/quality_control'
fastqc_res_path = quality_control_path + '/fastqc'
multiqc_res_path = quality_control_path + '/multiqc'
adaptor_removal_res_path = quality_control_path + '/adaptor'
trimming_res_path = quality_control_path + '/trimming'
src_path = rootDir + '/src'
file_handle_path = "${src_path}/file_handle/src/file_handle.py"
params.file_handle_version = "0.1.1"
params.name = "quality control analysis"
params.fastqc = "/usr/local/bin/fastqc"
params.fastqc_version = "0.11.5"
params.multiqc = "/usr/local/bin/multiqc"
params.multiqc_version = "1.0"
params.urqt = "/usr/local/bin/UrQt"
params.urqt_version = "d62c1f8"
params.cutadapt = "/usr/local/bin/cutadapt"
params.python2_version = "2.7"
params.python3_version = "3.5"
params.nextflow_version = "0.25.1"
params.gzip = "/usr/local/bin/gzip"
params.pigz = "/usr/local/bin/pigz"
params.pigz_version = "2.3.4"
params.do_adapter_removal = true

if (config.docker.enabled || params.global_executor == "sge") {
  file_handle_path = "/usr/bin/local/file_handle.py"
} else {
  if( !file(params.fastqc).exists() ) exit 1, "fastqc binary not found at: ${params.fastqc}"
  if( !file(params.multiqc).exists() ) exit 1, "multiqc binary not found at: ${params.multiqc}"
  if( !file(params.urqt).exists() ) exit 1, "UrQt binary not found at: ${params.urqt}"
  if( !file(params.cutadapt).exists() ) exit 1, "cutadapt binary not found at: ${params.cutadapt}"
}
params.adaptor_removal = "cutadapt"
params.adaptor_sequence = "-a AGATCGGAAGAG -g CTCTTCCGATCT"
params.paired = true
if(params.paired){
  params.adaptor_sequence = "-a AGATCGGAAGAG -g CTCTTCCGATCT -A AGATCGGAAGAG -G CTCTTCCGATCT"
}
params.trimmer = "UrQt"
params.quality_threshold = 20
if(params.paired != true && params.paired != false){
   exit 1, "Invalid paired option: ${params.paired}. Valid options: 'true' or 'false'"
}

process_header = ""
file_handle_module = ""
fastqc_module = ""
multiqc_module = ""
urqt_module = ""
cutadapt_module = ""
pigz_module = ""
python2_unload_mnodule = ""
python3_unload_mnodule = ""
if (params.global_executor == 'sge'){
  process_header = params.pbs_header
  file_handle_path = "file_handle.py"
  file_handle_module = "module load file_handle/${params.file_handle_version}"
  fastqc_module = "module load FastQC/${params.fastqc_version}"
  multiqc_module = "module load python/${params.python2_version}"
  urqt_module = "module load UrQt/${params.urqt_version}"
  cutadapt_module = "module load python/${params.python2_version}"
  pigz_module = "module load pigz/${params.pigz_version}"
  python2_unload_mnodule = "module unload python/${params.python2_version}"
  python3_unload_mnodule = "module unload python/${params.python3_version}"
}

log.info params.name
log.info "============================================"
log.info "executor: ${params.global_executor}"
log.info "fastq files : ${params.fastq_files}"
log.info "paired files : ${params.paired}"
if (params.paired) {
  log.info "file names are expected to end in the format *_{1,2}.fastq*."
  log.info "or *_R{1,2}.fastq*."
  log.info "otherwise the pairs will not be paired for the analysis"
}
log.info "fastqc : ${params.fastqc}"
log.info "multiqc : ${params.multiqc}"
log.info "results folder : ${results_path}"
log.info "trimmer : ${params.trimmer}"
if (params.trimmer == 'cutadapt') {
  log.info "cutadapt path : ${params.cutadapt}"
}else{
  log.info "UrQt path : ${params.urqt}"
}
gzip = ""
if(config.docker.enabled || params.global_executor == 'sge' || file(params.pigz).exists() ){
  gzip = params.pigz
  process get_pigz_version {
    echo true
    input:
      val params.pigz
    script:
"""
${pigz_module}
echo "\$(${params.pigz} --version)" &> grep pigz
"""
  }
  gzip = params.pigz
}else{
  log.info "pigz not found at ${params.pigz} using gzip"
  if( !file(params.gzip).exists() ) exit 1, "gzip binary not found at: ${params.gzip}"
  process get_gzip_version {
    echo true
    input:
      val params.gzip
    script:
"""
echo "\$(${params.gzip} --version)"
"""
  }
  gzip = params.gzip
}
log.info "gz software: ${gzip}"
log.info "number of cpu : ${params.cpu}"
log.info "\n"

fastq_files = Channel.fromFilePairs( params.fastq_files, size: -1)
  .ifEmpty { exit 1, "Cannot find any fastq files matching: ${params.fastq_files}" }

process get_file_name {
  tag "${tagname}"
  input:
    set val(fastq_name), file(reads) from fastq_files
  output:
    file "*.fastq.gz" into dated_file_names
  script:
    cmd_date = "${file_handle_path} -c -e -f"
    gzip_arg = ""
    if (gzip == params.pigz) { gzip_arg = "-p ${task.cpus}" }
    cmd_gzip = "${gzip} ${gzip_arg} -c "
    if (!(
      reads =~ /^.*\.fastq$/ || \
      reads =~ /^.*\.fastq\.gz$/
      )) {
      exit 1, "Can only work with fastq or fastq.gz files: ${reads}"
    }
    single = reads instanceof Path
    if (single) {
      tagname = (reads =~ /(.*\/){0,1}(.*)\.fastq(\.gz){0,1}/)[0][2]
      reads_0 = (reads =~ /(.*\/){0,1}(.*)/)[0][2]
      if (reads =~ /.*\.gz/) {
"""
${process_header}
${file_handle_module}
${cmd_date} *.fastq.gz
"""
      } else {
        reads_0 = reads_0 + ".gz"
"""
${process_header}
${pigz_module}
${cmd_gzip} ${reads} > ${reads_0}
${file_handle_module}
${cmd_date} ${reads_0}
"""
      }
    } else {
      tagname = (reads[0] =~ /(.*\/){0,1}(.*)_(R){0,1}[0,1]\.fastq(\.gz){0,1}/)[0][2]
      reads_1 = (reads[0] =~ /(.*\/){0,1}(.*)/)[0][2]
      reads_2 = (reads[1] =~ /(.*\/){0,1}(.*)/)[0][2]
      if (reads[0] =~ /.*\.gz/ || reads[1] =~ /.*\.gz/) {
"""
${process_header}
${file_handle_module}
${cmd_date} *.fastq.gz
"""
      } else {
        reads_1 = reads_1 + ".gz"
        reads_2 = reads_2 + ".gz"
"""
${process_header}
${pigz_module}
${cmd_gzip} ${reads[0]} > ${reads_1}
${cmd_gzip} ${reads[1]} > ${reads_2}
${file_handle_module}
${cmd_date} ${reads_1} ${reads_2}
"""
      }
    }
}

dated_file_names.into{
  fastqc_input;
  adaptor_rm_input
}

process fastqc {
  tag "${tagname}"
  publishDir "${fastqc_res_path}", mode: 'copy'
  input:
     file reads from fastqc_input
  output:
    file "*.{zip,html}" into fastqc_output
  script:
    cmd_date = "${file_handle_path} -c -f *.{zip,html}"
    single = reads instanceof Path
    if (single) {
      tagname = (reads =~ /(.*\/){0,1}(.*)\.fastq(\.gz){0,1}/)[0][2]
"""
${process_header}
${fastqc_module}
${params.fastqc} --quiet --outdir ./ ${reads}
${file_handle_module}
${cmd_date}
"""
    } else {
      tagname = (reads[0] =~ /(.*\/){0,1}(.*)_(R){0,1}[0,1]\.fastq(\.gz){0,1}/)[0][2]
"""
${process_header}
${fastqc_module}
${params.fastqc} --quiet --threads ${task.cpus} --outdir ./ ${reads[0]} ${reads[1]}
${file_handle_module}
${cmd_date}
"""
    }
}

process adaptor_removal {
  tag "${tagname}"
  publishDir "${adaptor_removal_res_path}", mode: 'copy'
  input:
    file reads from adaptor_rm_input
  output:
    file "*_cut*.fastq.gz" into adaptor_rm_output
    file "*_report.txt" into adaptor_rm_log
  script:
    cmd_date = "${file_handle_path} -c -f *_cut*.fastq.gz *_report.txt"
    single = reads instanceof Path
    if (single) {
      tagname = (reads =~ /(.*\/){0,1}(.*)\.fastq(\.gz){0,1}/)[0][2]
      reads_0 = tagname
      if (params.do_adapter_removal) {
"""
${process_header}
${cutadapt_module}
${params.cutadapt} ${params.adaptor_sequence} -o ${reads_0}_cut.fastq.gz ${reads} > ${tagname}_report.txt
${python2_unload_mnodule}
${file_handle_module}
${cmd_date}
"""
      } else {
"""
${process_header}
cp ${reads} ${reads_0}_cut.fastq.gz
printf 'No adapter removal done' > ${tagname}_report.txt
${python2_unload_mnodule}
${file_handle_module}
${cmd_date}
"""
      }
    } else {
      tagname = (reads[0] =~ /(.*\/){0,1}(.*)_(R){0,1}[0,1]\.fastq(\.gz){0,1}/)[0][2]
      reads_1 = "${tagname}_cut_R1.fastq.gz"
      reads_2 = "${tagname}_cut_R2.fastq.gz"
      if (params.do_adapter_removal) {
"""
${process_header}
${cutadapt_module}
${params.cutadapt} ${params.adaptor_sequence} -o ${reads_1} -p ${reads_2} ${reads[0]} ${reads[1]} > ${tagname}_report.txt
${python2_unload_mnodule}
${file_handle_module}
${cmd_date}
"""
      } else {
"""
${process_header}
cp ${reads[0]} ${reads_1}
cp ${reads[1]} ${reads_2}
printf 'No adapter removal done' > ${tagname}_report.txt
${python2_unload_mnodule}
${file_handle_module}
${cmd_date}
"""
      }
    }
}

process trimming {
  tag "${tagname}"
  publishDir "${trimming_res_path}", mode: 'copy'
  input:
    file reads from adaptor_rm_output
  output:
    file "*_trim*.fastq.gz" into trimming_output
    file "*_report.txt" into trimming_log
  script:
    cmd_date = "${file_handle_path} -r -f *_trim*.fastq.gz *_report.txt"
    single = reads instanceof Path
    if (single) {
      basename = (reads =~ /(.*)_cut\.fastq\.gz/)[0][1]
      tagname = basename
      if (params.trimmer == "cutadapt") {
"""
${process_header}
${cutadapt_module}
${params.cutadapt} -q ${params.quality_threshold},${params.quality_threshold} -o ${basename}_trim.fastq.gz ${reads} > ${basename}_cutadapt_report.txt
${python2_unload_mnodule}
${file_handle_module}
${cmd_date}
"""
      }else{
"""
${process_header}
${urqt_module}
${params.urqt} --m ${task.cpus} --t ${params.quality_threshold} --gz --in ${reads} --out ${basename}_trim.fastq.gz > ${basename}_UrQt_report.txt
${file_handle_module}
${cmd_date}
"""
      }
    } else {
      tagname = (reads[0] =~ /(.*)_cut_(R){0,1}[12]\.fastq\.gz/)[0][1]
      basename_1 = "${tagname}_trim_R1.fastq.gz"
      basename_2 = "${tagname}_trim_R2.fastq.gz"
      if (params.trimmer == "cutadapt") {
"""
${process_header}
${cutadapt_module}
${params.cutadapt} -q ${params.quality_threshold},${params.quality_threshold} -o ${basename_1} -p ${basename_2} ${reads[0]} ${reads[1]} > ${tagname}_cutadapt_report.txt
${python2_unload_mnodule}
${file_handle_module}
${cmd_date}
"""
      }else{
"""
${process_header}
${urqt_module}
${params.urqt} --m ${task.cpus} --t ${params.quality_threshold} --gz --in ${reads[0]} --inpair ${reads[1]} --out ${basename_1} --outpair ${basename_2} > ${tagname}_UrQt_report.txt
${file_handle_module}
${cmd_date}
"""
      }
    }
}

process fastqc_trimmed {
  tag "${tagname}"
  publishDir "${fastqc_res_path}", mode: 'copy'
  input:
    file reads from trimming_output
  output:
    file "*_fastqc.{zip,html}" into fastqc_trimmed_output
  script:
    cmd_date = "${file_handle_path} -c -f *.{zip,html}"
    single = reads instanceof Path
    if (single) {
      basename = (reads =~ /(.*)_trim\.fastq\.gz/)[0][1]
      tagname = basename
"""
${process_header}
${fastqc_module}
${params.fastqc} --quiet --outdir ./ ${reads}
${file_handle_module}
${cmd_date}
"""
    } else {
      tagname = (reads[0] =~ /(.*)_trim_(R){0,1}[12]\.fastq.gz/)[0][1]
"""
${process_header}
${fastqc_module}
${params.fastqc} --quiet --threads ${task.cpus} --outdir ./ ${reads[0]} ${reads[1]}
${file_handle_module}
${cmd_date}
"""
    }
}

process multiqc {
    publishDir "${multiqc_res_path}", mode: 'copy'
    input:
      file fast_qc_results from fastqc_output.collect()
      file fastqc_trimmed_results from fastqc_trimmed_output.collect()
    output:
      file "*multiqc_{report,data}*" into multiqc_report
    script:
"""
${multiqc_module}
${params.multiqc} -f .
${python2_unload_mnodule}
${file_handle_module}
${file_handle_path} -f multiqc_report.html -r
${file_handle_path} -f multiqc_data -r
"""
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
