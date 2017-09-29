#!/usr/bin/env nextflowsrc_path = rootDir + '/src'

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
rootDir = (baseDir =~ /(.*)src/)[0][1]
results_path = rootDir + 'results'
src_path = rootDir + '/src'


// test software path
class software_path {
  def params

  def call(params, docker, src_path) {
    try{
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
    try{
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
    try{
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
    try{
      return "${this.params.file_handle} -c -e -r -f ${file}"
    } catch (e) {
      println "error in software_path.cmd_date() ${e}"
    }
  }

  def cmd_unsalt_file(file) {
    try{
      file = (reads =~ /d(\d{4}_\d{2}_\d{2}_.*)/)[0][1]
    } catch (e) {
      println "error in software_path.cmd_unsalt_file() ${e}"
    }
  }

  def test_exist_fastq(file) {
    try{
      if (!(file ==~ /^.*\.fastq$/ || file ==~ /^.*\.fastq\.gz$/)) {
        exit 1, "Can only work with fastq or fastq.gz files: ${file}"
      }
    } catch (e) {
      println "error in software_path.test_exist_fastq() ${e}"
    }
  }

  def test_fastq(file) {
    try{
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
    try{
      if (this.test_single(file)) {
        return (file =~ /(.*\/){0,1}(.*)\.fastq(\.gz){0,1}/)[0][2]
      } else {
        return (file[0] =~ /(.*\/){0,1}(.*)_(R){0,1}[0,1]\.fastq(\.gz){0,1}/)[0][2]
      }
    } catch (e) {
      println "error in software_path.get_tagname() ${e}"
    }
  }

  def test_gz(file) {
    try{
      return file ==~ /^.*\.gz$/
    } catch (e) {
      println "error in software_path.test_gz() ${e}"
    }
  }

  def test_single(file) {
    try{
      return file instanceof Path
    } catch (e) {
      println "error in software_path.test_single() ${e}"
    }
  }
}

class modularity {
  def todo = [:]

  def call(path) {
    todo.list = split(path.params.list, '+')
  }
}

path = new software_path()
path(params, config.docker.enabled == true, src_path)
config.docker.runOptions = "--cpus=\"${path.params.cpu}\" --memory=\"${path.params.memory}\""

// load data
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
    template "${rootDir}src/func/get_file_name.sh"
}

//
//
// process fastqc {
//   tag "${tagname}"
//   publishDir "${fastqc_res_path}", mode: 'copy'
//   input:
//      file reads from fastqc_input
//   output:
//     file "*.{zip,html}" into fastqc_output
//   script:
//     cmd_date = "${file_handle_path} -c -f *.{zip,html}"
//     single = reads instanceof Path
//     if (single) {
//       tagname = (reads =~ /(.*\/){0,1}(.*)\.fastq(\.gz){0,1}/)[0][2]
// """
// ${process_header}
// ${fastqc_module}
// ${params.fastqc} --quiet --outdir ./ ${reads}
// ${file_handle_module}
// ${cmd_date}
// """
//     } else {
//       tagname = (reads[0] =~ /(.*\/){0,1}(.*)_(R){0,1}[0,1]\.fastq(\.gz){0,1}/)[0][2]
// """
// ${process_header}
// ${fastqc_module}
// ${params.fastqc} --quiet --threads ${task.cpus} --outdir ./ ${reads[0]} ${reads[1]}
// ${file_handle_module}
// ${cmd_date}
// """
//     }
// }
