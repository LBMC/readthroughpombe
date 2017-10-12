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

// test software path
class software_path {
  def params = [:]

  def call(params, docker, src_path) {
    try {
      for (param in params){
        this.params[param.key] = param.value
      }
      if (this.params.global_executor ==  'sge') {
          println 'executor : sge\n'
          this.params.gz = params.pigz
          def program_name
          for (param in params) {
            if (param.key ==~ /^.*_module$/) {
              program_name = (param.key =~ /^(.*)_module$/)[0][1]
              this.params[params.key] = this.params[params.key]+this.params["${program_name}_version"]
            }
          }
      } else {
        if (docker){
          println 'executor : docker\n'
          this.params.gz = params.pigz
        } else {
          println 'executor : local\n'
          println "path : ${src_path}"
          println "target : ${src_path}/file_handle/src/file_handle.py"
          println this.params.file_handle
          this.params.file_handle = "${src_path}/file_handle/src/file_handle.py"
          println this.params.file_handle
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
      if (file(this.params.pigz).exists()) {
        this.params.gz = this.params.pigz
      } else {
        this.params.gz = this.params.gzip
      }
    } catch (e) {
      println "error in software_path.test_pigz() ${e}"
    }
  }

  def cmd_gz(cpu, file) {
    try {
      def cmd = ""
      if (this.test_gz(file)) {
        return cmd
      } else {
        if (this.params.gz ==~ /.*pigz$/) {
          cmd = "${this.params.gz} -p ${cpu}"
        }else{
          cmd = "${this.params.gz}"
        }
        if (this.test_single(file)) {
          return "cat ${file} | ${cmd} -c > ${file}.gz"
        } else {
          return "cat ${file[0]} | ${cmd} -c > ${file[0]}.gz && \
          cat ${file[1]} | ${cmd} -c > ${file[1]}.gz"
        }
      }
    } catch (e) {
      println "error in software_path.cmd_gz() ${e}"
    }
  }

  def cmd_gz(cpu) {
    try {
      def cmd
      if (this.params.gz ==~ /.*pigz$/) {
        cmd = "${this.params.gz} -p ${cpu}"
      }else{
        cmd = "${this.params.gz}"
      }
      return "${cmd} -c"
    } catch (e) {
      println "error in software_path.cmd_gz() ${e}"
    }
  }

  def cmd_ungz(cpu, file) {
    try {
      def cmd
      file = this.unsalt_file_name(file)
      if (this.params.gz ==~ /.*pigz$/) {
        cmd = "${this.params.gz} -p ${cpu}"
      }else{
        cmd = "${this.params.gz}"
      }
      return "cat ${file} | ${cmd} -d -c"
    } catch (e) {
      println "error in software_path.cmd_gz() ${e}"
    }
  }

  def cmd_ungz(cpu) {
    try {
      def cmd
      if (this.params.gz ==~ /.*pigz$/) {
        cmd = "${this.params.gz} -p ${cpu}"
      }else{
        cmd = "${this.params.gz}"
      }
      return "${cmd} -d -c"
    } catch (e) {
      println "error in software_path.cmd_gz() ${e}"
    }
  }

  def cmd_date(file) {
    try {
      return "${this.params.file_handle} -e -r -f ${file}"
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
      println "error in software_path.test_unsalt() ${e}"
    }
  }

  def ungz_file_name(file){
    try {
      file = this.unsalt_file_name(file)
      if (this.test_single(file)) {
        file = (file =~ /^(.*)\.gz$/)[0][1]
      } else {
        file[0] = (file[0] =~ /^(.*)\.gz$/)[0][1]
        file[1] = (file[1] =~ /^(.*)\.gz$/)[0][1]
      }
      return file
    } catch (e) {
      println "error in software_path.ungz_file_name() ${e}"
    }
  }

  def unsalt_file_name(file){
    try {
      if (this.test_single(file)) {
        file = (file =~ /^d{0,1}(\d{4}_\d{2}_\d{2}_.*)/)[0][1]
      } else {
        file[0] = (file[0] =~ /^d{0,1}(\d{4}_\d{2}_\d{2}_.*)/)[0][1]
        file[1] = (file[1] =~ /^d{0,1}(\d{4}_\d{2}_\d{2}_.*)/)[0][1]
      }
      return file
    } catch (e) {
      println "error in software_path.unsalt_file_name() ${e}"
    }
  }

  def cmd_unsalt_file(file) {
    try {
      if (this.test_unsalt(file)) {
        return """
find . -name "d*" | \
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
      def cmd = "${this.params.cutadapt}"
      def tagname = this.get_tagname(file)
      if (this.test_single(file)) {
        return "${cmd} ${this.params.adaptor_sequence_single} -o ${tagname}_cut.fastq.gz ${file} > ${tagname}_report.txt"
      } else {
        return "${cmd} ${this.params.adaptor_sequence_paired} -o ${tagname}_cut_R1.fastq.gz -p ${tagname}_cut_R2.fastq.gz ${file[0]} ${file[1]} > ${tagname}_report.txt"
      }
    } catch (e) {
      println "error in software_path.cmd_adaptor_removal() ${e}"
    }
  }

  def cmd_cutadapt(cpu, file){
    try {
      file = this.unsalt_file_name(file)
      def tagname = this.get_tagname(file)
      if (this.test_single(file)) {
        return "${this.params.cutadapt} -q ${this.params.quality_threshold},${this.params.quality_threshold} -o ${tagname}_trim.fastq.gz ${file} > ${tagname}_cutadapt_report.txt"
      } else {
        return "${this.params.cutadapt} -q ${this.params.quality_threshold},${this.params.quality_threshold} -o ${tagname}_trim_R1.fastq.gz -p ${tagname}_trim_R2.fastq.gz ${file[0]} ${file[1]} > ${tagname}_cutadapt_report.txt"
      }
    } catch (e) {
      println "error in software_path.cmd_cutadapt() ${e}"
    }
  }


  def cmd_urqt(cpu, file){
    try {
      file = this.unsalt_file_name(file)
      def tagname = this.get_tagname(file)
      if (this.test_single(file)) {
        return "${this.params.urqt} --m ${cpu} --t ${this.params.quality_threshold} --gz --in ${file} --out ${tagname}_trim.fastq.gz > ${tagname}_UrQt_report.txt"
      } else {
        return "${this.params.urqt} --m ${cpu} --t ${this.params.quality_threshold} --gz --in ${file[0]} --inpair ${file[1]} --out ${tagname}_trim_R1.fastq.gz --outpair ${tagname}_trim_R2.fastq.gz > ${tagname}_UrQt_report.txt"
      }
    } catch (e) {
      println "error in software_path.cmd_urqt() ${e}"
    }
  }

  def basename_index(file){
    try {
      if ( file in nextflow.util.BlankSeparatedList ) {
        return (file[0] =~ /^(.*\.index).*$/)[0][1]
      } else {
        return (file =~ /^(.*\.index).*$/)[0][1]
      }
    } catch (e) {
      println "error in software_path.basename_index() ${e}"
    }
  }

  def cmd_kallisto(cpu, index, fastq){
    try {
      fastq = this.unsalt_file_name(fastq)
      index = this.unsalt_file_name(index)
      def tagname = this.get_tagname(fastq)
      def basename_index = this.basename_index(index)
      if (this.test_single(fastq)) {
        return "${this.params.kallisto} quant -i ${basename_index} -t ${cpu} --single ${this.params.kallisto_parameters} -l ${this.params.mean} -s ${this.params.sd} -o ./ ${fastq} &> ${tagname}_kallisto_report.txt"
      } else {
        return "${this.params.kallisto} quant -i ${basename_index} -t ${cpu} ${this.params.kallisto_parameters} -o ./ ${fastq[0]} ${fastq[1]} &> ${tagname}_kallisto_report.txt"
      }
    } catch (e) {
      println "error in software_path.cmd_kallisto() ${e}"
    }
  }

  def cmd_bowtie2(cpu, index, fastq, todo){
    try {
      fastq = this.unsalt_file_name(fastq)
      index = this.unsalt_file_name(index)
      def tagname = this.get_tagname(fastq)
      def basename_index = this.basename_index(index)
      if (this.test_single(fastq)) {
        return "${this.params.bowtie2} ${this.params.bowtie2_parameters} -p ${cpu} -x ${basename_index} -U ${fastq} 2> ${tagname}_bowtie2_report.txt | samtools view -Sb - > ${tagname}.bam"
      } else {
        return "${this.params.bowtie2} ${this.params.bowtie2_parameters} -p ${cpu} -x ${basename_index} -1 ${fastq[0]} -2 ${fastq[1]} 2> ${tagname}_bowtie2_report.txt | samtools view -Sb - > ${tagname}.bam"
      }
    } catch (e) {
      println "error in software_path.cmd_bowtie2() ${e}"
    }
  }

  def cmd_rsem_bowtie2_index(cpu, reference, annotation){
    try {
      annotation = this.unsalt_file_name(annotation)
      reference = this.unsalt_file_name(reference)
      def tagname = this.get_tagname(reference)
      def bowtie2_path = '$(which ' + this.params.bowtie2 + ' | sed \'s/bowtie2$//g\')'
      def cmd = "${this.params.rsem}-prepare-reference -p ${cpu} --bowtie2 --bowtie2-path ${bowtie2_path}"
      if (this.params.rsem_parameters_indexing != "") {
        cmd = "${cmd} ${this.params.rsem_parameters_indexing} "
      }
      def cmd_annotation = "--gff3 ${annotation}"
      if(annotation ==~ /.*\.gtf$/){
        cmd_annotation = "--gtf ${annotation}"
      }
      return "${cmd} ${cmd_annotation} ${reference} ${tagname}.index > ${tagname}_rsem_bowtie2_report.txt"
    } catch (e) {
      println "error in software_path.cmd_rsem_bowtie2_index() ${e}"
    }
  }

  def cmd_rsem_bowtie2(cpu, fastq, index){
    try {
      fastq = this.unsalt_file_name(fastq)
      index = this.unsalt_file_name(index)
      def tagname = this.get_tagname(fastq)
      def basename_index = this.basename_index(index)
      def bowtie2_path = '$(which ' + this.params.bowtie2 + ' | sed \'s/bowtie2$//g\')'
      if (this.test_single(fastq)) {
        return "${this.params.rsem}-calculate-expression --bowtie2 --bowtie2-path ${bowtie2_path} ${this.params.rsem_parameters_quantif} -p ${cpu} ${fastq} ${basename_index} ${tagname} > ${tagname}_rsem_bowtie2_report.txt"
      } else {
        return "${this.params.rsem}-calculate-expression --bowtie2 --bowtie2-path ${bowtie2_path} --bowtie2-sensitivity-level \"very_sensitive\" ${this.params.rsem_parameters_quantif} -p ${cpu} --paired-end ${fastq[0]} ${fastq[1]} ${basename_index} ${tagname} > ${tagname}_rsem_bowtie2_report.txt"
      }
    } catch (e) {
      println "error in software_path.cmd_rsem_bowtie2() ${e}"
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

  def test_exist_bam(file) {
    try {
      if (!(file ==~ /^.*\.bam$/)) {
        exit 1, "Can only work with fastq or fastq.gz files: ${file}"
      }
    } catch (e) {
      println "error in software_path.test_exist_bam() ${e}"
    }
  }

  def test_bam(file) {
    try {
      if (this.test_single(file)) {
        this.test_exist_bam(file)
      } else {
        for (bam in file){
          this.test_exist_bam(bam)
        }
      }
    } catch (e) {
      println "error in software_path.test_fastq() ${e}"
    }
  }

  def test_exist_fasta(file) {
    try {
      if (!(file ==~ /^.*\.fasta$/ || file ==~ /^.*\.fasta\.gz$/)) {
        exit 1, "Can only work with fasta or fasta.gz files: ${file}"
      }
    } catch (e) {
      println "error in software_path.test_exist_fasta() ${e}"
    }
  }

  def test_fasta(file) {
    try {
      if (this.test_single(file)) {
        this.test_exist_fasta(file)
      } else {
        for (fasta in file){
          this.test_exist_fasta(fasta)
        }
      }
    } catch (e) {
      println "error in software_path.test_fasta() ${e}"
    }
  }

  def test_exist_annot(file) {
    try {
      if (!(file ==~ /^.*\.gff(.gz){0,1}$/ || \
        file ==~ /^.*\.gtf(.gz){0,1}$/ || \
        file ==~ /^.*\.gff3(.gz){0,1}$/)) {
        exit 1, "Can only work with gff(3) or gtf files: ${file}"
      }
    } catch (e) {
      println "error in software_path.test_exist_annot() ${e}"
    }
  }

  def test_annot(file) {
    try {
      if (this.test_single(file)) {
        this.test_exist_annot(file)
      } else {
        return false
      }
    } catch (e) {
      println "error in software_path.test_annot() ${e}"
    }
  }

  def get_tagname(file){
    try {
      if (this.test_single(file)) {
        if (file ==~ /^.*\.fastq$/ || file ==~ /^.*\.fastq\.gz$/) {
          return (file =~ /(.*\/){0,1}d{0,1}(.*)\.fast[aq](\.gz){0,1}/)[0][2]
        }
        if (file ==~ /^.*\.bam$/) {
          return (file =~ /(.*\/){0,1}d{0,1}(.*)\.bam/)[0][2]
        }
        if (file ==~ /^.*\.fasta$/ || file ==~ /^.*\.fasta\.gz$/) {
          return (file =~ /(.*\/){0,1}d{0,1}(.*)\.fasta(\.gz){0,1}/)[0][2]
        }
      } else {
        return (file[0] =~ /(.*\/){0,1}d{0,1}(.*)_(R){0,1}[0,1]\.fastq(\.gz){0,1}/)[0][2]
      }
      return file
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
      if ( file in nextflow.util.BlankSeparatedList ) {
        return false
      } else {
        return true
      }
    } catch (e) {
      println "error in software_path.test_single() ${e}"
    }
  }
}

class modularity {
  def todo = [
    'fastqc_raw' : ['fastqc'],
    'adaptor_rm' : ['cutadapt'],
    'trimming' : ['urqt', 'cutadapt'],
    'fastqc_trim' : ['fastqc'],
    'multiqc_qc' : ['multiqc'],
    'split_ref' : ['split_ref'],
    'indexing' : ['none'],
    'mapping' : ['bowtie2', 'kallisto'],
    'quantification' : ['htseq', 'rsem'],
    'multiqc_mapping' : ['multiqc']
  ]
  def path = List

  def call(path) {
    this.path = path
    def todo_list = path.params.todo.replaceAll("\\s","").tokenize('+')
    def job_number = 0
    for (job in this.todo) {
      if (todo_list[job_number] in this.todo[job.key]) {
        this.todo[job.key] = todo_list[job_number]
        job_number += 1
      } else {
        this.todo[job.key] = 'none'
      }
    }
    this.do_split_ref()
    this.do_indexing_ref()
    this.display_task()
  }

  def do_split_ref(){
    if (this.reference() && this.annotation()) {
      if (this.todo['mapping'] in ['kallisto']) {
        this.todo['split_ref'] = 'split_fasta'
        this.todo['quantification'] = 'kallisto'
      }
      if (this.todo['mapping'] in ['bowtie2'] && this.todo['quantification'] in ['rsem']) {
        this.todo['split_ref'] = 'none'
        this.todo['mapping'] = 'rsem_bowtie2'
        this.todo['quantification'] = 'rsem_bowtie2'
      }
    }
  }

  def do_indexing_ref(){
    if (this.reference() && (this.mapping() || this.pseudomapping())){
      this.todo['indexing'] = this.todo['mapping']
    } else {
      this.todo['indexing'] = 'none'
    }
  }

  def display_task(){
    def job_number = 0
    for (job in this.todo) {
      if (this.todo[job.key] != "none") {
        println "${job.key} : ${this.todo[job.key]}"
        job_number += 1
      }
    }
    if (job_number > 0) {
      println "${job_number} tasks to do."
    } else {
      println "error: no task found for ${todo_list}"
    }
  }

  def fastqc_raw(){
    return this.todo['fastqc_raw'] != 'none'
  }
  def adaptor_removal(){
    return this.todo['adaptor_rm'] != 'none'
  }
  def trimming(){
    return this.todo['trimming'] != 'none'
  }
  def fastqc_trim(){
    return this.todo['fastqc_trim'] != 'none'
  }
  def multiqc_qc(){
    return this.todo['multiqc_qc'] != 'none'
  }
  def reference(){
    return this.path.params.fasta != ""
  }
  def annotation(){
    return this.path.params.annot != ""
  }
  def split_ref(){
    return this.todo['split_ref'] != 'none'
  }
  def indexing(){
    if (this.reference()) {
      return this.todo['indexing'] != 'none'
    }
    return false
  }
  def load_index(){
    return this.path.params.index != ""
  }
  def mapping(){
    if (this.reference()) {
      return !(this.todo['mapping'] in ['none', 'kallisto', 'rsem_bowtie2'])
    }
    return false
  }
  def pseudomapping(){
    if (this.reference()) {
      return this.todo['mapping'] in ['kallisto', 'rsem_bowtie2']
    }
    return false
  }
  def quantification(){
    if (this.reference()) {
      return !(this.todo['quantification'] in ['none', 'kallisto', 'rsem_bowtie2'])
    }
    println "warning: no fasta provided skipping quantification"
    return false
  }
  def multiqc_mapping(){
    if (this.reference) {
      return this.todo['multiqc_mapping'] != 'none'
    }
    return false
  }
}

// global variables
root_path = (baseDir =~ /(.*)src/)[0][1]
results_path = root_path + 'results'
src_path = root_path + 'src'
path = new software_path()
path(params, config.docker.enabled == true, src_path)
todo = new modularity()
todo(path)
config.docker.runOptions = "--cpus=\"${path.params.cpu}\" --memory=\"${path.params.memory}\""

fastqc_output = Channel.create()
////////////////////////////////// load fastq //////////////////////////////////
if(path.params.fastq == ""){exit 1, "missing params \"--fastq\""}
log.info "fastq files : ${path.params.fastq}"
fastq_files = Channel.fromFilePairs(path.params.fastq, size: -1)
  .ifEmpty { exit 1, "Cannot find any fastq files matching: ${path.params.fastq}" }

process get_fastq_name {
  tag "${tagname}"
  echo path.params.verbose
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

////////////////////////////////// load fasta //////////////////////////////////
log.info "reference files : ${path.params.fasta}"
if (todo.reference()) {
  fasta_files = Channel.fromPath( path.params.fasta )
  process get_fasta_name {
    tag "${tagname}"
    echo path.params.verbose
    input:
      file reference from fasta_files
    output:
      file "d*.fasta.gz" into dated_fasta_files
    script:
      path.test_fasta(reference)
      tagname = reference
      file = reference
      template "${src_path}/func/get_file_name.sh"
  }
}

////////////////////////////////// load annot //////////////////////////////////
log.info "annotation files : ${path.params.annot}"
if (todo.annotation()) {
  annot_files = Channel.fromPath( path.params.annot )
  process get_fasta_name {
    tag "${tagname}"
    echo path.params.verbose
    input:
      file annotation from annot_files
    output:
      file "d*" into dated_annot_files
    script:
      path.test_annot(annotation)
      tagname = annotation
      file = annotation
      template "${src_path}/func/get_file_name.sh"
  }
}

////////////////////////////////// load index //////////////////////////////////
if (todo.load_index()) {
  log.info "index files : ${path.params.index}"
  index_files = Channel.fromPath( path.params.index )
  process get_fasta_name {
    tag "${tagname}"
    echo path.params.verbose
    input:
      file index from index_files
    output:
      file "d*" into index_file_1
    script:
      path.test_annot(index)
      tagname = index
      file = index
      template "${src_path}/func/get_file_name_nogz.sh"
  }
}

/////////////////////////// fastqc on raw fastq ////////////////////////////////
if (todo.fastqc_raw()) {
  dated_fastq_files.into{
    fastqc_raw_input;
    fastq_file_1
  }
  process fastqc_raw {
    tag "${tagname}"
    echo path.params.verbose
    publishDir "${results_path}/quality_control/fastqc", mode: 'copy'
    input:
       file reads from fastqc_raw_input
    output:
      file "*.{zip,html}" into fastqc_raw_output
    script:
      path.test_fastq(reads)
      tagname = path.get_tagname(reads)
      file = reads
      template "${src_path}/func/quality_control/fastqc.sh"
  }
} else {
  dated_fastq_files.into{
    fastqc_raw_output;
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
    echo path.params.verbose
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
    echo path.params.verbose
    publishDir "${results_path}/quality_control/trimming", mode: 'copy'
    input:
      file reads from trimming_input
    output:
      file "*_trim*.fastq.gz" into fastq_file_3
      file "*_report.txt" into trimming_log
    script:
        path.test_fastq(reads)
        tagname = path.get_tagname(reads)
        file = reads
        template "${src_path}/func/quality_control/${todo.todo.trimming}.sh"
  }
} else {
  fastq_file_2.set{
    fastq_file_3
  }
}

/////////////////////////// fastqc on trim fastq ////////////////////////////////
if (todo.fastqc_trim()) {
  fastq_file_3.into{
    fastqc_trim_input;
    fastq_file_4
  }
  process fastqc_trim {
    tag "${tagname}"
    echo path.params.verbose
    publishDir "${results_path}/quality_control/fastqc", mode: 'copy'
    input:
       file reads from fastqc_trim_input
    output:
      file "*.{zip,html}" into fastqc_trim_output
    script:
      path.test_fastq(reads)
      tagname = path.get_tagname(reads)
      file = reads
      template "${src_path}/func/quality_control/fastqc.sh"
  }
} else {
  fastq_file_3.into{
    fastqc_trim_output;
    fastq_file_4
  }
}

//////////////////////////////// multiqc on QC /////////////////////////////////
if (todo.multiqc_qc()) {
  fastqc_raw_output.into{
    multiqc_qc_input_raw;
    multiqc_mapping_input_raw
  }
  fastqc_trim_output.into{
    multiqc_qc_input_trim;
    multiqc_mapping_input_trim
  }
  process multiqc_qc {
    echo path.params.verbose
    publishDir "${results_path}/quality_control/multiqc_qc", mode: 'copy'
    input:
       file file_fastqc_raw from multiqc_qc_input_raw.collect()
       file file_fastqc_trim from multiqc_qc_input_trim.collect()
    output:
      file "*multiqc_*" into multiqc_report
    script:
      template "${src_path}/func/quality_control/multiqc.sh"
  }
} else {
  fastqc_raw_output.set{
    multiqc_mapping_input_raw
  }
  fastqc_trim_output.set{
    multiqc_mapping_input_trim
  }
}

///////////////////////////////// split ref ////////////////////////////////////
if (todo.split_ref()) {
  dated_annot_files.into{
    split_ref_annot;
    annot_file_1
  }
  process split_ref {
    tag "${tagname}"
    echo path.params.verbose
    publishDir "${results_path}/mapping/split_fasta", mode: 'copy'
    input:
       file reference from dated_fasta_files
       file annotation from split_ref_annot
    output:
      file "*_split.fasta.gz" into fasta_file_1
    script:
      path.test_fasta(reference)
      path.test_annot(annotation)
      tagname = path.get_tagname(reference)
      template "${src_path}/func/mapping/split_fasta.sh"
  }
} else {
  dated_fasta_files.set{
    fasta_file_1
  }
  dated_annot_files.set{
    annot_file_1
  }
}

////////////////////////////////// indexing ////////////////////////////////////
if (todo.indexing()) {
  annot_file_1.into {
    indexing_annot_file;
    annot_file_2
  }
  process indexing {
    tag "${tagname}"
    echo path.params.verbose
    publishDir "${results_path}/mapping/indexing", mode: 'copy'
    input:
       file reference from fasta_file_1
       file annotation from indexing_annot_file
    output:
      file "*index*" into index_file_1
      file "*_report.txt" into index_report
    script:
      path.test_fasta(reference)
      tagname = path.get_tagname(reference)
      template "${src_path}/func/mapping/${todo.todo.indexing}_index.sh"
  }
}else{
  annot_file_1.set {
    annot_file_2
  }
}

/////////////////////////////// pseudomapping //////////////////////////////////
if (todo.pseudomapping()) {
  process pseudomapping {
    tag "${tagname}"
    echo path.params.verbose
    publishDir "${results_path}/mapping/quantification", mode: 'copy'
    input:
       file index from index_file_1.first()
       file reads from fastq_file_4
    output:
      file "*" into quantification_file_1
    script:
      path.test_fastq(reads)
      tagname = path.get_tagname(reads)
      template "${src_path}/func/mapping/${todo.todo.mapping}.sh"
  }
}

/////////////////////////////////// mapping ////////////////////////////////////
if (todo.mapping()) {
  process mapping {
    tag "${tagname}"
    echo path.params.verbose
    publishDir "${results_path}/mapping/mapping", mode: 'copy'
    input:
      file index from index_file_1.first()
      file reads from fastq_file_4
    output:
      file "*.bam" into mapping_file_1
      file "*_report.txt" into mapping_report
    script:
      path.test_fastq(reads)
      tagname = path.get_tagname(reads)
      template "${src_path}/func/mapping/${todo.todo.mapping}.sh"
  }

  process sort_bam {
    tag "${tagname}"
    echo path.params.verbose
    publishDir "${results_path}/mapping/mapping", mode: 'copy'
    input:
      file bam from mapping_file_1
    output:
      file "*_sort.bam" into mapping_file_2
    script:
      path.test_bam(bam)
      tagname = path.get_tagname(bam)
      template "${src_path}/func/mapping/sort_bam.sh"
  }
}

/////////////////////////////// quantification /////////////////////////////////
if (todo.quantification()) {
  process quantification {
    tag "${tagname}"
    echo path.params.verbose
    publishDir "${results_path}/mapping/quantification", mode: 'copy'
    input:
       file bam from mapping_file_2
       file annotation from annot_file_2.first()
    output:
      file "*" into quantification_file_1
    script:
      path.test_bam(bam)
      path.test_annot(annotation)
      tagname = path.get_tagname(bam)
      template "${src_path}/func/mapping/${todo.todo.quantification}.sh"
  }
}
