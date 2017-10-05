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

// global variables
println baseDir
root_path = (baseDir =~ /(.*)src/)[0][1]
results_path = root_path + 'results'
src_path = root_path + '/src'


// test software path
class software_path {
  def params

  def call(params, docker, src_path) {
    try {
      this.params =  params
      if (this.params.global_executor ==  'sge') {
          println 'executor : sge'
          this.params.gz = params.pigz
      } else {
        if (docker){
          println 'executor : docker\n'
          this.params.gz = params.pigz
        } else {
          println 'executor : local\n'
          this.params.file_handle_path = "${src_path}/file_handle/src/file_handle.py"
          this.test_pigz()
        }
      }
      this.params.process_header = params.pbs_header
    } catch (e) {
      println "error in software_path.call() ${e}"
    }
  }

  def test_pigz() {
    try {
      if (file(params.pigz).exists()) {
        this.params.gz = params.pigz
      } else {
        this.params.gz = params.gzip
      }
    } catch (e) {
      println "error in software_path.test_pigz() ${e}"
    }
  }

  def cmd_gz(cpu, file) {
    try {
      def cmd
      if (this.params.gz ==~ /.*pigz$/) {
        cmd = "${this.params.gz} -p ${cpu}"
      }else{
        cmd = "${this.params.gz}"
      }
      if (this.test_single(file)) {
        return "cat ${file} | ${cmd} -c > ${file}.gz"
      } else {
        return "cat ${file[1]} | ${cmd} -c > ${file[0]}.gz && \
        cat ${file[1]} | ${cmd} -c > ${file[1]}.gz"
      }
    } catch (e) {
      println "error in software_path.cmd_gz() ${e}"
    }
  }

  def cmd_date(file) {
    try {
      return "${this.params.file_handle} -c -e -r -f ${file}"
    } catch (e) {
      println "error in software_path.cmd_date() ${e}"
    }
  }

  def test_unsalt(file){
    try {
      if (this.test_single(file)) {
        return file ==~ /^d\d{4}_\d{2}_\d{2}_.*/
      } else {
        return file[0] ==~ /^d\d{4}_\d{2}_\d{2}_.*/ ||
          file[1] ==~ /^d\d{4}_\d{2}_\d{2}_.*/
      }
    } catch (e) {
      println "error in software_path.cmd_unsalt_file() ${e}"
    }
  }

  def unsalt_file_name(file){
    try {
      if (this.test_single(file)) {
        file = (file =~ /^d(\d{4}_\d{2}_\d{2}_.*)/)[0][1]
      } else {
        file[0] = (file[0] =~ /^d(\d{4}_\d{2}_\d{2}_.*)/)[0][1]
        file[1] = (file[1] =~ /^d(\d{4}_\d{2}_\d{2}_.*)/)[0][1]
      }
      return file
    } catch (e) {
      println "error in software_path.cmd_unsalt_file() ${e}"
    }
  }

  def cmd_unsalt_file(file) {
    try {
      if (this.test_unsalt(file)) {
        return """
find . -type f -name "d*" | \
sed 's/^.\\/d//g' | \
awk '{system("mv d"\$0" "\$0)}'
"""
      } else {
        return ""
      }
    } catch (e) {
      println "error in software_path.cmd_unsalt_file() ${e}"
    }
  }

  def cmd_fastqc(cpu, file){
    try {
      file = this.unsalt_file_name(file)
      def cmd = "${this.params.fastqc} --quiet --threads ${cpu} --outdir ./"
      if (this.test_single(file)) {
        return "${cmd} ${file}"
      } else {
        return "${cmd} ${file[0]} ${file[1]}"
      }
    } catch (e) {
      println "error in software_path.cmd_fastqc() ${e}"
    }
  }

  def cmd_adaptor_removal(cpu, file){ try {
      file = this.unsalt_file_name(file)
      def cmd = "${params.cutadapt}"
      def tagname = this.get_tagname(file)
      if (this.test_single(file)) {
        return "${cmd} ${params.adaptor_sequence_single} -o ${tagname}_cut.fastq.gz ${file} > ${tagname}_report.txt"
      } else {
        return "${cmd} ${params.adaptor_sequence_paired} -o ${tagname}_cut_R1.fastq.gz -p ${tagname}_cut_R2.fastq.gz ${file[0]} ${file[1]} > ${tagname}_report.txt"
      }
    } catch (e) {
      println "error in software_path.cmd_fastqc() ${e}"
    }
  }

  def trimming_module(){
    try {
      switch(this.params.trimmer) {
        case 'cutadapt':
          return this.params.cutadapt_module
        break
        default:
          return this.params.urqt_module
        break
      }
    } catch (e) {
      println "error in software_path.trimming_module() ${e}"
    }
  }

  def cmd_trimming(cpu, file){
    try {
      file = this.unsalt_file_name(file)
      def tagname = this.get_tagname(file)
      if (this.test_single(file)) {
        switch(this.params.trimmer) {
          case 'cutadapt':
            return "${this.params.cutadapt} -q ${this.params.quality_threshold},${this.params.quality_threshold} -o ${basename}_trim.fastq.gz ${file} > ${basename}_cutadapt_report.txt"
          break
          default:
            return "${this.params.urqt} --m ${cpu} --t ${this.params.quality_threshold} --gz --in ${file} --out ${basename}_trim.fastq.gz > ${basename}_UrQt_report.txt"
          break
        }
      } else {
        switch(this.params.trimmer) {
          case 'cutadapt':
            return "${this.params.cutadapt} -q ${this.params.quality_threshold},${this.params.quality_threshold} -o ${tagname}_trim_R1.fastq.gz -p ${tagname}_trim_R2.fastq.gz ${file[0]} ${file[1]} > ${tagname}_cutadapt_report.txt"
          break
          default:
            return "${this.params.urqt} --m ${cpu} --t ${this.params.quality_threshold} --gz --in ${file[0]} --inpair ${file[1]} --out ${tagname}_trim_R1.fastq.gz --outpair ${tagname}_trim_R2.fastq.gz > ${tagname}_UrQt_report.txt"
          break
        }
      }
    } catch (e) {
      println "error in software_path.cmd_trimming() ${e}"
    }
  }

  def test_exist_fastq(file) {
    try {
      if (!(file ==~ /^.*\.fastq$/ || file ==~ /^.*\.fastq\.gz$/)) {
        exit 1, "Can only work with fastq or fastq.gz files: ${file}"
      }
    } catch (e) {
      println "error in software_path.test_exist_fastq() ${e}"
    }
  }

  def test_fastq(file) {
    try {
      if (this.test_single(file)) {
        this.test_exist_fastq(file)
      } else {
        this.test_exist_fastq(file[0])
        this.test_exist_fastq(file[1])
      }
    } catch (e) {
      println "error in software_path.test_fastq() ${e}"
    }
  }

  def get_tagname(file){
    try {
      if (this.test_single(file)) {
        return (file =~ /(.*\/){0,1}d{0,1}(.*)\.fastq(\.gz){0,1}/)[0][2]
      } else {
        return (file[0] =~ /(.*\/){0,1}d{0,1}(.*)_(R){0,1}[0,1]\.fastq(\.gz){0,1}/)[0][2]
      }
    } catch (e) {
      println "error in software_path.get_tagname() ${e}"
    }
  }

  def test_gz(file) {
    try {
      return file ==~ /^.*\.gz$/
    } catch (e) {
      println "error in software_path.test_gz() ${e}"
    }
  }

  def test_single(file) {
    try {
      return file instanceof Path
    } catch (e) {
      println "error in software_path.test_single() ${e}"
    }
  }
}

class modularity {
  def todo = [
    'fastqc_raw' : [true, 'fastqc'],
    'adaptor_rm' : [true, 'cutadapt'],
    'trimming' : [true, 'urqt'],
    'fastqc_trim' : [true, 'fastqc'],
    'multiqc_trim' : [true, 'multiqc'],
    'mapping' : [true, 'bowtie2'],
    'quantifying' : [true, 'htseq'],
    'multiqc_mapping' : [true, 'multiqc']
  ]

  def call(path) {
    this.todo['adaptor_rm'][1] = path.params.adaptor_removal
    this.todo['trimming'][1] = path.params.trimmer
    this.todo['mapping'][1] = path.params.mapper
    this.todo['quantifying'][1] = path.params.quantifier

    def todo_list = path.params.todo.replaceAll("\\s","").tokenize('+')
    def job_number = 0
    for (job in this.todo) {
      if (!todo_list[job_number] == job.key[1]) {
        this.todo[job.key][0] = false
      } else {
        job_number += 1
      }
    }
  }

  def fastqc_raw(){
    return this.todo['fastqc_raw'][0]
  }
  def adaptor_removal(){
    return this.todo['adaptor_rm'][0]
  }
  def trimming(){
    return this.todo['trimming'][0]
  }
  def fastqc_trim(){
    return this.todo['fastqc_trim'][0]
  }
  def multiqc_trim(){
    return this.todo['multiqc_trim'][0]
  }
  def indexing(){
    return this.todo['indexing'][0]
  }
  def mapping(){
    return this.todo['mapping'][0]
  }
  def quantifying(){
    return this.todo['quantifying'][0]
  }
  def multiqc_mapping(){
    return this.todo['multiqc_mapping'][0]
  }
}

path = new software_path()
path(params, config.docker.enabled == true, src_path)
todo = new modularity()
todo(path)
config.docker.runOptions = "--cpus=\"${path.params.cpu}\" --memory=\"${path.params.memory}\""

/////////////////////////////// load fastq /////////////////////////////////////
fastq_files = Channel.fromFilePairs( params.fastq, size: -1)
  .ifEmpty { exit 1, "Cannot find any fastq files matching: ${params.fastq}" }

process get_fastq_name {
  tag "${tagname}"
  input:
    set val(fastq_name), file(reads) from fastq_files
  output:
    file "d*.fastq.gz" into dated_fastq_files
  script:
    path.test_fastq(reads)
    tagname = path.get_tagname(reads)
    file = reads
    template "${src_path}/func/get_file_name.sh"
}

/////////////////////////// fastqc on raw fastq ////////////////////////////////
if (todo.fastqc_raw()) {
  dated_fastq_files.into{
    fastqc_raw_input;
    fastq_file_1
  }
  process fastqc_raw {
    tag "${tagname}"
    echo true
    publishDir "${results_path}/quality_control/fastqc", mode: 'copy'
    input:
       file reads from fastqc_raw_input
    output:
      file "*.{zip,html}" into fastqc_output
    script:
      path.test_fastq(reads)
      tagname = path.get_tagname(reads)
      file = reads
      template "${src_path}/func/quality_control/fastqc.sh"
  }
} else {
  dated_fastq_files.set{
    fastq_file_1
  }
}

//////////////////////////////// adaptor_removal////////////////////////////////
if (todo.adaptor_removal()) {
  fastq_file_1.set{
    adaptor_rm_input;
  }
  process adaptor_removal {
    tag "${tagname}"
    echo true
    publishDir "${results_path}/quality_control/adaptor", mode: 'copy'
    input:
      file reads from adaptor_rm_input
    output:
      file "*_cut*.fastq.gz" into fastq_file_2
      file "*_report.txt" into adaptor_rm_log
    script:
        path.test_fastq(reads)
        tagname = path.get_tagname(reads)
        file = reads
        template "${src_path}/func/quality_control/adaptor_removal.sh"
  }
} else {
  fastq_file_1.set{
    fastq_file_2
  }
}

//////////////////////////////// trimming //////////////////////////////////////
if (todo.trimming()) {
  fastq_file_2.set{
    trimming_input;
  }
  process trimming {
    tag "${tagname}"
    publishDir "${results_path}/quality_control/trimming", mode: 'copy'
    input:
      file reads from trimming_input
    output:
      file "*_trimming*.fastq.gz" into fastq_file_3
      file "*_report.txt" into trimming_log
    script:
        path.test_fastq(reads)
        tagname = path.get_tagname(reads)
        file = reads
        template "${src_path}/func/quality_control/trimming.sh"
  }
} else {
  fastq_file_2.set{
    fastq_file_3
  }
}
