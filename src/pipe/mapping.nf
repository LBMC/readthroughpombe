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
params.reference = ""
params.fastq_files = ""

log.info params.name
log.info "============================================"
if(params.fastq_files == ""){exit 1, "missing params \"--fastq_files\""}
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
if(params.reference == ""){exit 1, "missing params \"--reference\""}
log.info "reference files : ${params.reference}"
reference_names = Channel.fromPath( params.reference)
  .ifEmpty { exit 1, "Cannot find any reference file matching: ${params.reference}" }
if(params.annotation == ""){exit 1, "missing params \"--annotation\""}
log.info "annotation files : ${params.annotation}"
annotation_names = Channel.fromPath(params.annotation)
  .ifEmpty { exit 1, "Cannot find any annotation file matching: ${params.annotation}" }
switch(params.mapper) {
  case "salmon":
    log.info "salmon path : ${params.salmon}"
    if( !file(params.salmon).exists() ) exit 1, "salmon binary not found at: ${params.salmon}"
    process get_salmon_version {
      echo true
      input:
        val params.salmon
      script:
      """
      ${params.salmon} --version &> salmon_version.txt
      echo "salmon \$(cat salmon_version.txt)"
      """
    }
  break
  case "kallisto":
    log.info "kallisto path : ${params.kallisto}"
    if( !file(params.kallisto).exists() ) exit 1, "kallisto binary not found at: ${params.kallisto}"
    process get_kallisto_version {
      echo true
      input:
        val params.kallisto
      script:
      """
      echo "\$(${params.kallisto} version)"
      """
    }
  break
  case "bowtie2":
    log.info "bowtie2 path : ${params.bowtie2}"
    if( !file(params.bowtie2).exists() ) exit 1, "bowtie2 binary not found at: ${params.bowtie2}"
    if( !file(params.bowtie2+"-build").exists() ) exit 1, "bowtie2-build binary not found at: ${params.bowtie2}-build"
    process get_bowtie2_version {
      echo true
      input:
        val params.bowtie2
      script:
      """
      echo "\$(${params.bowtie2} --version)"
      """
    }
  break
  default:
  exit 1, "Invalid paired option: ${params.mapper}. Valid options: 'salmon', 'kallisto' or 'bowtie2'"
  break
}
log.info "bedtools path : ${params.bedtools}"
if( !file(params.bedtools).exists() ) exit 1, "bedtools binary not found at: ${params.bedtools}"
log.info "samtools path : ${params.samtools}"
process get_bedtools_version {
  echo true
  input:
    val params.bedtools
  script:
  """
  echo "\$(${params.bedtools} --version)"
  """
}
if( !file(params.samtools).exists() ) exit 1, "samtools binary not found at: ${params.samtools}"
log.info "results folder : ${results_path}"
process get_samtools_version {
  echo true
  input:
    val params.samtools
  script:
  """
  echo "\$(${params.samtools} --version)"
  """
}
log.info "\n"

process get_file_name_fastq {
  tag "${tagname}"
  input:
    val fastq_name from fastq_names
  output:
    file "*${fastq_name_root}" into dated_fastq_names
  when:
    if (params.paired) {
      (fastq_name[1][0] =~ /^.*\.fastq$/ || fastq_name[1][0] =~ /^.*\.fastq\.gz$/) && (fastq_name[1][1] =~ /^.*\.fastq$/ || fastq_name[1][1] =~ /^.*\.fastq\.gz$/)
    } else {
      fastq_name =~ /^.*\.fastq$/ || fastq_name =~ /^.*\.fastq\.gz$/
    }
  script:
  if (params.paired) {
    tagname = (fastq_name[1][0].baseName =~ /^(.*)_(R){0,1}[0,1]\.fastq(\.gz){0,1}/)[0][1]
    fastq_name_root = fastq_name[0] + "*"
    """
    ${src_path}/func/file_handle.py -f ${fastq_name[1][0]} ${fastq_name[1][1]} -c -e | \
    awk '{system("ln -s "\$0" ."); print(\$0)}'
    """
  } else {
    tagname = file(fastq_name).baseName
    fastq_name_root = fastq_name.name
    """
    ${src_path}/func/file_handle.py -f ${fastq_name} -c -e | \
    awk '{system("ln -s "\$0" ."); print(\$0)}'
    """
  }
}

process get_file_name_reference {
  tag "${tagname}"
  input:
    val reference_name from reference_names
  output:
    file "*${reference_name_root}" into dated_reference_names
  when:
    reference_name =~ /^.*\.fasta$/ || reference_name =~ /^.*\.fasta\.gz$/ || reference_name =~ /^.*\.fasta\.index$/ || reference_name =~ /^.*\.fasta\.gz\.index$/
  script:
    tagname = file(reference_name).baseName
    reference_name_root = reference_name.name
    """
    ${src_path}/func/file_handle.py -f ${reference_name} -c -e | \
    awk '{system("ln -s "\$0" ."); print(\$0)}'
    """
}

dated_reference_names
  .into{ dated_reference_names_split; dated_reference_names_mapping }

process get_file_name_annotation {
  tag "${tagname}"
  input:
    val annotation_name from annotation_names
  output:
    file "*${file_name_root}" into dated_annotation_names
  when:
    params.annotation != "" & (annotation_name =~ /^.*\.gtf$/ || annotation_name =~ /^.*\.bed$/ || annotation_name =~ /^.*\.gff$/ || annotation_name =~ /^.*\.vcf$/)
  script:
    tagname = file(annotation_name).baseName
    file_name_root = annotation_name.name
    """
    ${src_path}/func/file_handle.py -f ${annotation_name} -c -e | \
    awk '{system("ln -s "\$0" ."); print(\$0)}'
    """
}

dated_annotation_names
  .into{ dated_annotation_names_split; dated_annotation_names_quantification }

if(params.mapper in ["salmon", "kallisto"]){
  process split_ref {
    tag "${tagname}"
    publishDir "${index_res_path}", mode: 'copy'
    cpu = 12
    input:
      file reference_name from dated_reference_names_split
      file annotation_name from dated_annotation_names_split
    output:
      file "*_split.fasta" into indexing_input
    when:
      params.mapper in ["salmon", "kallisto"]
    script:
    if ( reference_name ==~ /.*\.index/){
      exit 1, "Cannot split an index file with a annotation file. Provide a fasta file instead of  ${params.reference}"
    }
    basename = (reference_name =~ /(.*)(\.gff){0,1}(\.bed){0,1}(\.vcf){0,1}(\.gtf){0,1}/)[0][1]
    basename_fasta = (reference_name =~ /(.*)\.fasta/)[0][1]
    tagname = basename
    """
    cat ${basename} | gunzip > ${basename_fasta}.fasta
    ${params.bedtools} getfasta -fi ${basename_fasta}.fasta -bed ${annotation_name} -fo ${basename_fasta}_split.fasta
    ${src_path}/func/file_handle.py -f *_split.fasta -r
    ls -lh
    """
  }
}else{
  dated_reference_names_split.into{ indexing_input }
}

process indexing {
  tag "${tagname}"
  publishDir "${index_res_path}", mode: 'copy'
  cpu = 12
  input:
    file index_name from indexing_input
  output:
    file "*.index*" into indexing_output
    file "*_report.txt*" into indexing_log
  script:
  basename = (index_name =~ /(.*\.fasta)(\.index){0,1}/)[0][1]
  tagname = basename
  if ( index_name ==~ /.*\.index/) {
    log.info "index file found. Skipping indexing step"
  } else {
    switch(params.mapper) {
      case "kallisto":
        """
        ${params.kallisto} index -k 31 --make-unique -i ${basename}.index ${index_name} > ${basename}_kallisto_report.txt
        ${src_path}/func/file_handle.py -f *.index -r
        """
      break
      case "bowtie2":
        """
        ${params.bowtie2}-build --threads ${task.cpu} ${index_name} ${basename}.index > ${basename}_bowtie2_report.txt
        ${src_path}/func/file_handle.py -f *.index* -r
        """
      default:
        """
        ${params.salmon} index -p ${task.cpu} -t ${index_name} -i ${basename}.index --type quasi -k 31 &> ${basename}_salmon_report.txt
        ${src_path}/func/file_handle.py -f *.index* -r
        """
      break
    }
  }
}

process mapping {
  tag "${tagname}"
  if(params.mapper in ["salmon", "kallisto"]){
    publishDir "${counts_res_path}", mode: 'copy'
  }else{
    publishDir "${bams_res_path}", mode: 'copy'
  }
  cpu = 12
  input:
    file index_name from indexing_output
    file fastq_name from dated_fastq_names
  output:
    if(!params.mapper in ["salmon", "kallisto"]){
      file "*.bam" into mapping_output
    }
    file "*_report.txt" into mapping_log
  script:
  if (params.paired) {
    name = (fastq_name[0] =~ /(.*)\_(R){0,1}[12]\.fastq(\.gz){0,1}/)[0][1]
    tagname = name
    basename_1 = (fastq_name[0] =~ /(.*)\.fastq(\.gz){0,1}/)[0][1]
    basename_2 = (fastq_name[1] =~ /(.*)\.fastq(\.gz){0,1}/)[0][1]
    switch(params.mapper) {
      case "kallisto":
        """
        ${params.kallisto} quant -i ${index_name} -t ${task.cpu} ${params.kallisto_parameters} -o ./ ${fastq_name[0]} ${fastq_name[1]} &> ${name}_kallisto_report.txt
        mv abundance.tsv ${name}.counts
        mv run_info.json ${name}_info.json
        mv abundance.h5 ${name}.h5
        ${src_path}/func/file_handle.py -f *_report.txt *.counts *.json *.h5 -r
        """
      break
      case "bowtie2":
        """
        ${params.bowtie2} ${params.bowtie2_parameters} -p ${task.cpu} -x *.index* -1 ${fastq_name[0]} -2 ${fastq_name[1]} 2> ${name}_bowtie2_report.txt | samtools view -Sb - > ${name}.bam
        ${src_path}/func/file_handle.py -f * -r
        """
      break
      default:
        """
        ${params.salmon} quant -i ${index_name} -p ${task.cpu} ${params.salmon_parameters} -1 ${fastq_name[0]} -2 ${fastq_name[1]} -o ${name}.counts > ${name}_salmon_report.txt
        ${src_path}/func/file_handle.py -f * -r
        """
      break
    }
  } else {
    name = (fastq_name =~ /(.*)\.fastq(\.gz){0,1}/)[0][1]
    tagname = name
    switch(params.mapper) {
      case "kallisto":
        kallisto_parameters = params.kallisto_parameters + " -l ${params.mean} -s ${params.sd}"
        """
        ${params.kallisto} quant -i ${index_name} -t ${task.cpu} --single ${kallisto_parameters} -o ./ ${fastq_name} > ${name}_report.txt
        mv abundance.tsv ${name}.counts
        mv run_info.json ${name}_info.json
        mv abundance.h5 ${name}.h5
        ${src_path}/func/file_handle.py -f * -r
        """
      break
      case "bowtie2":
        """
        ${params.bowtie2} ${params.bowtie2_parameters} -p ${task.cpu} -x ${index_name} -U ${fastq_name} 2> ${name}_report.txt | samtools view -Sb - > ${name}.bam
        ${src_path}/func/file_handle.py -f * -r
        """
      break
      default:
        salmon_parameters = params.salmon_parameters + " --fldMean ${params.mean} --fldSD ${params.sd}"
        """
        ${params.salmon} quant -i ${index_name} -p ${task.cpu} ${salmon_parameters} -r ${fastq_name} -o ${name}.counts > ${name}_report.txt
        ${src_path}/func/file_handle.py -f * -r
        """
      break
    }
  }
}

mapping_log.println()
