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
rootDir = (baseDir =~ /(.*)src\/pipe/)[0][1]
results_path = rootDir + 'results'
src_path = rootDir + '/src'

// test software path
class software_path {
  def process_header = ""
  def file_handle_path = ""
  def file_handle_module = ""
  def fastqc_module = ""
  def multiqc_module = ""
  def urqt_module = ""
  def cutadapt_module = ""
  def pigz_module = ""
  def python2_unload_module = ""
  def python3_unload_module = ""

  def call(params) {
    try{
      switch(params.global_executor) {
        case 'sge':
          println 'executor : sge'
          process_header = params.pbs_header
          this.file_handle_path = "file_handle.py"
          this.file_handle_module = "module load file_handle/${params.file_handle_version}"
          this.fastqc_module = "module load FastQC/${params.fastqc_version}"
          this.multiqc_module = "module load python/${params.python2_version}"
          this.urqt_module = "module load UrQt/${params.urqt_version}"
          this.cutadapt_module = "module load python/${params.python2_version}"
          this.pigz_module = "module load pigz/${params.pigz_version}"
          this.python2_unload_module = "module unload python/${params.python2_version}"
          this.python3_unload_module = "module unload python/${params.python3_version}"
        break
        case 'docker':
          println 'executor : docker\n'
          this.file_handle_path = "/usr/bin/local/file_handle.py"
        break
        default:
          println 'executor : local\n'
          this.file_handle_path = "${src_path}/file_handle/src/file_handle.py"
        break
      }
    } catch (e) {
      println "error in software_path.call() ${e}"
    }
  }
}

path = new software_path()
path(params)

// template "${rootDir}src/func/evaluate_parameters.groovy"
