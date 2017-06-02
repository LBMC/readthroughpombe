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
bams_res_path = quality_control_path + '/bams'
index_res_path = quality_control_path + '/bams'
counts_res_path = quality_control_path + '/counts'
src_path = rootDir + '/src'

params.name = "mapping and quantification analysis"
params.mapper = "salmon"
params.paired = true
params.salmon = "/usr/bin/salmon"
params.salmon_parameters = "--useVBOpt --numBootstraps 100 --seqBias --gcBias --posBias --libType MU"
params.kallisto = "/usr/local/bin/kallisto"
params.kallisto_parameters = "--bias --bootstrap-samples 100"
params.bowtie2 = "/usr/bin/bowtie2"
params.bowtie2_parameters = "--very-sensitive"
params.bedtools = "/usr/bin/bedtools"
params.samtools = "/usr/bin/samtools"
params.mean = 200
params.sd = 20
params.annotation = ""

log.info params.name
log.info "============================================"
log.info "fastq files : ${params.fastq_files}"
log.info "paired files : ${params.paired}"
if(params.paired != true && params.paired != false){
   exit 1, "Invalid paired option: ${params.paired}. Valid options: 'true' or 'false'"
}
if (params.paired) {
  log.info "file names are expected to end in the format *_{1,2}.fastq*."
  log.info "or *_R{1,2}.fastq*."
  log.info "otherwise the pairs will not be paired for the analysis"
  fastq_names = Channel.fromFilePairs( params.fastq_files, size:2)
    .ifEmpty { exit 1, "Cannot find any fastq files matching: ${params.fastq_files}" }
}else{
  log.info "mean fragment length : ${params.mean}"
  log.info "standar deviation fragment length : ${params.sd}"
  fastq_names = Channel.fromPath( params.fastq_files )
    .ifEmpty { exit 1, "Cannot find any fastq files matching: ${params.fastq_files}" }
}
log.info "reference files : ${params.reference}"
reference_names = Channel.fromPath( params.reference)
  .ifEmpty { exit 1, "Cannot find any reference file matching: ${params.reference}" }
if(params.annotation != ""){
  log.info "annotation files : ${params.annotation}"
  annotation_name = Channel.fromPath(params.annotation)
    .ifEmpty { exit 1, "Cannot find any annotation file matching: ${params.annotation}" }
}
switch(params.mapper) {
  case "salmon":
    log.info "salmon path : ${params.salmon}"
    if( !file(params.salmon).exists() ) exit 1, "salmon binary not found at: ${params.salmon}"
  break
  case "kallisto":
    log.info "kallisto path : ${params.kallisto}"
    if( !file(params.kallisto).exists() ) exit 1, "kallisto binary not found at: ${params.kallisto}"
  break
  case "bowtie2":
    log.info "bowtie2 path : ${params.bowtie2}"
    if( !file(params.bowtie2).exists() ) exit 1, "bowtie2 binary not found at: ${params.bowtie2}"
    if( !file(params.bowtie2+"-build").exists() ) exit 1, "bowtie2-build binary not found at: ${params.bowtie2}-build"
  break
  default:
  exit 1, "Invalid paired option: ${params.mapper}. Valid options: 'salmon' or 'kallisto'"
  break
}
log.info "bedtools path : ${params.bedtools}"
if( !file(params.bedtools).exists() ) exit 1, "bedtools binary not found at: ${params.bedtools}"
log.info "samtools path : ${params.samtools}"
if( !file(params.samtools).exists() ) exit 1, "samtools binary not found at: ${params.samtools}"
log.info "results folder : ${results_path}"
log.info "\n"

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

process get_file_name_reference {
  tag "${tagname}"
  input:
    val file_name from reference_names
  output:
    file "*${file_name_root}" into dated_reference_names
  when:
    file_name =~ /^.*\.fasta$/ || file_name =~ /^.*\.fasta\.gz$/ || file_name =~ /^.*\.fasta\.index$/ || file_name =~ /^.*\.fasta\.gz\.index$/
  script:
    tagname = file(file_name).baseName
    file_name_root = file_name.name
    """
    ${src_path}/func/file_handle.py -f ${file_name} -c -e | \
    awk '{system("ln -s "\$0" ."); print(\$0)}'
    """
}

process get_file_name_annotation {
  tag "${tagname}"
  input:
    val file_name from annotation_names
  output:
    file "*${file_name_root}" into dated_ann ${params.bowtie2_parameters}otation_names
  when:
    params.annotation != "" & (file_name =~ /^.*\.gtf$/ || file_name =~ /^.*\.bed$/ || file_name =~ /^.*\.gff$/ || file_name =~ /^.*\.vcf$/)
  script:
    tagname = file(file_name).baseName
    file_name_root = file_name.name
    """
    ${src_path}/func/file_handle.py -f ${file_name} -c -e | \
    awk '{system("ln -s "\$0" ."); print(\$0)}'
    """
}

process split_ref {
  tag "${tagname}"
  publishDir "${index_res_path}", mode: 'copy'
  cpu = 12
  input:
    file file_name from dated_reference_names
    file annotation_name from dated_annotation_names
  output:
    file "*.fasta" into split_dated_reference_names
  when:
    params.mapper == "salmon" || params.mapper == "kallisto"
  script:
  if ( file_name ==~ /.*\.index/){
    exit 1, "Cannot split an index file with a annotation file. Provide a fasta file instead of  ${params.reference}"
  }
  basename = (file_name =~ /(.*(\.gff){0,1}(\.bed){0,1}(\.vcf){0,1}(\.gtf){0,1}/)[0][1]
  basename_fasta = (file_name =~ /(.*\.fasta)/)[0][1]
  tagname = basename
  """
  ${params.bedtools} getfasta -fi ${file_name} -bed ${annotation_name} -fo ${basename_fasta}_split.fasta
  ${src_path}/func/file_handle.py -f *.fasta -r
  """
}

process indexing {
  tag "${tagname}"
  publishDir "${index_res_path}", mode: 'copy'
  cpu = 12
  input:
    if(params.mapper == "salmon" || params.mapper == "kallisto"){
      file file_name from split_dated_reference_names
    }else{
      file file_name from dated_reference_names
    }
  output:
    if(params.mapper == "salmon" || params.mapper == "kallisto"){
      file "*.index" into indexing_output
    }else{
      file "*${base_name}.index*" into indexing_output
    }
  script:
  basename = (file_name =~ /(.*\.fasta)(\.index){0,1}/)[0][1]
  tagname = basename
  if ( file_name ==~ /.*\.index/) {
    log.info "index file found. Skipping indexing step"
  } else {
    switch(params.mapper) {
      case "kallisto":
        """
        ${params.kallisto} index -k 31 -i ${file_name}.index ${file_name}
        ${src_path}/func/file_handle.py -f *.index -r
        """
      break
      case "bowtie2":
        """
        ${params.bowtie2}-build --threads ${task.cpu} ${file_name} ${base_name}.index
        ${src_path}/func/file_handle.py -f ${base_name}.index* -r
        """
      default:
        """
        ${params.salmon} index -p ${task.cpu} -t ${file_name} -i ${file_name}.index --type quasi -k 31
        ${src_path}/func/file_handle.py -f *.index -r
        """
      break
    }
  }
}

process mapping {
  tag "${tagname}"
  if(params.mapper == "salmon" || params.mapper == "kallisto"){
    publishDir "${counts_res_path}", mode: 'copy'
  }else{
    publishDir "${bams_res_path}", mode: 'copy'
  }
  cpu = 12
  input:
    file index_name from indexing_output
    file file_name from dated_fastq_names
  output:
    if(params.mapper == "salmon" || params.mapper == "kallisto"){
    }else{
      file "*.bam" into mapping_output
    }
    file "*_report.txt" into mapping_log
  script:
  if (params.paired) {
    name = (file_name[0] =~ /(.*)\_(R){0,1}[12]\.fastq(\.gz){0,1}/)[0][1]
    tagname = name
    basename_1 = (file_name[0] =~ /(.*)\.fastq(\.gz){0,1}/)[0][1]
    basename_2 = (file_name[1] =~ /(.*)\.fastq(\.gz){0,1}/)[0][1]
    switch(params.mapper) {
      case "kallisto":
        """
        ${params.kallisto} quant -i ${index_name} -t ${task.cpu} ${params.kallisto_parameters} -o ./ ${file_name[0]} ${file_name[1]} > ${name}_report.txt
        mv abundance.tsv ${name}.counts
        mv run_info.json ${name}_info.json
        mv abundance.h5 ${name}.h5
        ${src_path}/func/file_handle.py -f * -r
        """
      break
      case "bowtie2":
        """
        ${params.bowtie2} ${params.bowtie2_parameters} -p ${task.cpu} -x *.index* -1 ${file_name[0]} -2 ${file_name[1]} | samtools view -Sb - > ${name}.bam
        """
      break
      default:
        """
        ${params.salmon} quant -i ${index_name} -p ${task.cpu} ${params.salmon_parameters} -1 ${file_name[0]} -2 ${file_name[1]} -o ${name}.counts > ${name}_report.txt
        ${src_path}/func/file_handle.py -f * -r
        """
      break
    }
  } else {
    name = (file_name =~ /(.*)\.fastq(\.gz){0,1}/)[0][1]
    tagname = name
    switch(params.mapper) {
      case "kallisto":
        kallisto_parameters = params.kallisto_parameters + " -l ${params.mean} -s ${params.sd}"
        """
        ${params.kallisto} quant -i ${index_name} -t ${task.cpu} --single ${kallisto_parameters} -o ./ ${file_name} > ${name}_report.txt
        mv abundance.tsv ${name}.counts
        mv run_info.json ${name}_info.json
        mv abundance.h5 ${name}.h5
        ${src_path}/func/file_handle.py -f * -r
        """
      break
      case "bowtie2":
        """
        ${params.bowtie2} ${params.bowtie2_parameters} -p ${task.cpu} -x ${index_name} -U ${file_name} | samtools view -Sb - > ${name}.bam
        """
      break
      default:
        salmon_parameters = params.salmon_parameters + " --fldMean ${params.mean} --fldSD ${params.sd}"
        """
        ${params.salmon} quant -i ${index_name} -p ${task.cpu} ${salmon_parameters} -r ${file_name} -o ${name}.counts > ${name}_report.txt
        ${src_path}/func/file_handle.py -f * -r
        """
      break
    }
  }
}
