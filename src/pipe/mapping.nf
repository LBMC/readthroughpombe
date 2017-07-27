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
index_res_path = quality_control_path + '/index'
counts_res_path = quality_control_path + '/counts'
src_path = rootDir + '/src'
file_handle_path = "${src_path}/file_handle/src/file_handle.py"

params.name = "mapping and quantification analysis"
params.mapper = "salmon"
params.paired = true
params.salmon = "/usr/local/bin/salmon"
params.salmon_version = "0.8.2"
params.salmon_parameters = "--useVBOpt --numBootstraps 100 --seqBias --gcBias --posBias --libType IU"
params.kallisto = "/usr/local/bin/kallisto"
params.kallisto_version = "0.43.1"
params.kallisto_parameters = "--bias --bootstrap-samples 100"
params.bowtie2 = "/usr/local/bin/bowtie2"
params.bowtie2_version = "2.3.2"
params.bowtie2_parameters = "--very-sensitive"
params.bedtools = "/usr/bin/bedtools"
params.bedtools_version = "2.25.0"
params.samtools = "/usr/local/bin/samtools"
params.samtools_version = "1.5"
params.quantifier = "htseq"
params.htseq = "/usr/local/bin/htseq-count"
params.htseq_version = "0.8.0"
params.htseq_parameters = "--mode=intersection-nonempty -a 10 -s no -t exon -i gene_id"
params.rsem = "/usr/local/bin/rsem"
params.rsem_version =
params.rsem_parameters = "1.3.0"
params.gzip = "/usr/bin/gzip"
params.pigz = "/usr/bin/pigz"
params.pigz_version = "2.3.4"
params.mean = 200
params.sd = 20
params.annotation = ""
params.reference = ""
params.fastq_files = ""
params.cpu = 12
params.file_handle_version = "0.1.1"
params.python2_version = "2.7"
params.python3_version = "3.5"
params.nextflow_version = "0.25.1"

if (config.docker.enabled || params.global_executor == "sge") {
  file_handle_path = "/usr/bin/local/file_handle.py"
}

process_header = ""
file_handle_module = ""
pigz_module = ""
salmon_module = ""
kallisto_module = ""
bowtie2_module = ""
bedtools_module = ""
samtools_module = ""
htseq_module = ""
rsem_module = ""
python2_unload_mnodule = ""
python3_unload_mnodule = ""
if (params.global_executor == 'sge'){
  process_header = params.pbs_header
  file_handle_path = "file_handle.py"
  file_handle_module = "module load file_handle/${params.file_handle_version}"
  pigz_module = "module load pigz/${params.pigz_version}"
  salmon_module = "module load Salmon/${params.salmon_version}"
  kallisto_module = "module load Kallisto/${params.kallisto_version}"
  bowtie2_module = "module load Bowtie2/${params.bowtie2_version}"
  bedtools_module = "module load BEDtools/${params.bedtools_version}"
  samtools_module = "module load SAMtools/${params.samtools_version}"
  htseq_module = "module load HTSeq/${params.htseq_version}"
  rsem_module = "module load RSEM/${params.rsem_version}"
  python2_unload_mnodule = "module unload python/${params.python2_version}"
  python3_unload_mnodule = "module unload python/${params.python3_version}"
}

log.info params.name
log.info "============================================"
if(params.fastq_files == ""){exit 1, "missing params \"--fastq_files\""}
log.info "fastq files : ${params.fastq_files}"
fastq_files = Channel.fromFilePairs( params.fastq_files, size: -1 )
  .ifEmpty { exit 1, "Cannot find any fastq files matching: ${params.fastq_files}" }

if(params.reference == ""){exit 1, "missing params \"--reference\""}
log.info "reference files : ${params.reference}"
reference_names = Channel.fromPath( params.reference )
  .ifEmpty { exit 1, "Cannot find any reference file matching: ${params.reference}" }

if(params.annotation == ""){exit 1, "missing params \"--annotation\""}
log.info "annotation files : ${params.annotation}"
annotation_names = Channel.fromPath( params.annotation )
  .ifEmpty { exit 1, "Cannot find any annotation file matching: ${params.annotation}" }

log.info "mapper : ${params.mapper}"
switch(params.mapper) {
  case "salmon":
    log.info "salmon path : ${params.salmon}"
    if( !config.docker.enabled && !params.global_executor == "sge" && !file(params.salmon).exists() ) exit 1, "salmon binary not found at: ${params.salmon}"
    process get_salmon_version {
      echo true
      input:
        val params.salmon
      script:
      """
      ${salmon_module}
      ${params.salmon} --version &> salmon_version.txt
      echo "salmon \$(cat salmon_version.txt)"
      """
    }
  break
  case "kallisto":
    log.info "kallisto path : ${params.kallisto}"
    if( !config.docker.enabled && !params.global_executor == "sge" && !file(params.kallisto).exists() ) exit 1, "kallisto binary not found at: ${params.kallisto}"
    process get_kallisto_version {
      echo true
      input:
        val params.kallisto
      script:
      """
      ${kallisto_module}
      echo "\$(${params.kallisto} version)"
      """
    }
  break
  case "bowtie2":
    log.info "bowtie2 path : ${params.bowtie2}"
    if( !config.docker.enabled && !params.global_executor == "sge" && !file(params.bowtie2).exists() ) exit 1, "bowtie2 binary not found at: ${params.bowtie2}"
    if( !config.docker.enabled && !params.global_executor == "sge" && !file(params.bowtie2+"-build").exists() ) exit 1, "bowtie2-build binary not found at: ${params.bowtie2}-build"
    process get_bowtie2_version {
      echo true
      input:
        val params.bowtie2
      script:
      """
      ${bowtie2_module}
      echo "\$(${params.bowtie2} --version | grep "bowtie")"
      """
    }
  break
  default:
  exit 1, "Invalid mapper option: ${params.mapper}. Valid options: 'salmon', 'kallisto' or 'bowtie2'"
  break
}
switch(params.quantifier) {
  case "rsem":
    log.info "rsem path : ${params.rsem}"
    if( !config.docker.enabled && !params.global_executor == "sge" && !file(params.rsem+"-prepare-reference").exists() ) exit 1, "rsem binaries not found at: ${params.rsem}-prepare-reference"
    if( !config.docker.enabled && !params.global_executor == "sge" && !file(params.rsem+"-calculate-expression").exists() ) exit 1, "rsem binaries not found at: ${params.rsem}-calculate-expression"
    if(!params.mapper in ["bowtie2"]) {
      exit 1, "RSEM can only work with bowtie2"
    }
    process get_rsem_version {
      echo true
      input:
        val params.rsem
      script:
      """
      ${rsem_module}
      echo "\$(${params.rsem}-calculate-expression --version)"
      """
    }
  break
  case "htseq":
    log.info "htseq path : ${params.htseq}"
    if( !config.docker.enabled && !params.global_executor == "sge" && !file(params.htseq).exists() ) exit 1, "htseq binaries not found at: ${params.htseq}"
    process get_htseq_version {
      echo true
      input:
        val params.htseq
      script:
      """
      ${htseq_module}
      echo "\$(${params.htseq} -h | grep "version")"
      """
    }
  break
  default:
  exit 1, "Invalid quantifier option: ${params.quantifier}. Valid options: 'rsem' or 'htseq'"
  break
}
log.info "bedtools path : ${params.bedtools}"
if( !config.docker.enabled && !params.global_executor == "sge" && !file(params.bedtools).exists() ) exit 1, "bedtools binary not found at: ${params.bedtools}"
log.info "samtools path : ${params.samtools}"
process get_bedtools_version {
  echo true
  input:
    val params.bedtools
  script:
  """
  ${bedtools_module}
  echo "\$(${params.bedtools} --version)"
  """
}
if( !config.docker.enabled && !params.global_executor == "sge" && !file(params.samtools).exists() ) exit 1, "samtools binary not found at: ${params.samtools}"
log.info "results folder : ${results_path}"
process get_samtools_version {
  echo true
  input:
    val params.samtools
  script:
  """
  ${samtools_module}
  echo "\$(${params.samtools} --version | grep "samtools")"
  """
}
gzip = ""
if( config.docker.enabled || params.global_executor == "sge" || file(params.pigz).exists() ){
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
mapper = params.mapper
rsem_parameters = params.rsem_parameters
if (params.mapper == "bowtie2" && params.quantifier == "rsem") {
  mapper = "bowtie2+rsem"
  bowtie2_path = (params.bowtie2 =~ /(.*)bowtie2/)[0][1]
  rsem_parameters = rsem_parameters + " --bowtie2 --bowtie2-path ${bowtie2_path}"
}

process get_fastq_name {
  tag "${tagname}"
  input:
    set val(fastq_name), file(reads) from fastq_files
  output:
    file "*.fastq.gz" into dated_fastq_names
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
        cp ${reads} ./${reads_0}
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
        cp ${reads[0]} ./${reads_1}
        cp ${reads[1]} ./${reads_2}
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

process get_file_name_reference {
  tag "${tagname}"
  input:
    file refs from reference_names
  output:
    file "*.fasta.gz" into dated_reference_names
  script:
    if (!(
      refs =~ /^.*\.fasta$/ || \
      refs =~ /^.*\.fasta\.gz$/
      )) {
      exit 1, "Can only work with fasta or fasta.gz files: ${refs}"
    }
    tagname = (refs =~ /(.*\/){0,1}(.*)\.fasta(\.gz){0,1}/)[0][2]
    reference = (refs =~ /(.*\/){0,1}(.*)/)[0][2]
    cmd_date = "${file_handle_path} -c -e -f"

    if (refs =~ /.*\.gz/) {
      """
      ${process_header}
      ${cmd_date} ${reference}
      """
    } else {
      reference = reference + ".gz"
      gzip_arg = ""
      if(gzip == params.pigz){gzip_arg = "-p ${task.cpus}"}
      """
      ${process_header}
      ${pigz_module}
      ${gzip} ${gzip_arg} -c ${refs} > ${reference}
      ${file_handle_module}
      ${cmd_date} ${reference}
      """
    }
}

dated_reference_names
  .into{
    dated_reference_names_split;
    dated_reference_names_mapping
  }

process get_file_name_annotation {
  tag "${tagname}"
  input:
    file annot from annotation_names
  output:
    file "*.{gtf,bed,gff,vcf}" into dated_annotation_names
  script:
    if (!(
      params.annotation != "" && \
        (annot =~ /^.*\.gtf$/ || \
          annot =~ /^.*\.bed$/ || \
          annot =~ /^.*\.gff$/ || \
          annot =~ /^.*\.vcf$/)
      )) {
      exit 1, "Can only work with gtf, bed, gff or vcf files: ${annot}"
    }
    tagname = (annot =~ /(.*\/){0,1}(.*)\.*/)[0][2]
    annotation = (annot =~ /(.*\/){0,1}(.*)/)[0][2]
    """
    ${file_handle_module}
    ${file_handle_path} -c -e -f ${annotation}
    """
}

dated_annotation_names
  .into{
    dated_annotation_names_split;
    dated_annotation_names_indexing;
    dated_annotation_names_quantification
  }

if(mapper in ["salmon", "kallisto"]){
  process split_ref {
    tag "${tagname}"
    publishDir "${index_res_path}", mode: 'copy'

    input:
      file reference_name from dated_reference_names_split
      file annotation_name from dated_annotation_names_split
    output:
      file "*_split.fasta" into indexing_input
    when:
      mapper in ["salmon", "kallisto"]
    script:
      if ( reference_name ==~ /.*\.index/) {
        exit 1, "Cannot split an index file with a annotation file. Provide a fasta file instead of  ${params.reference}"
      }
      basename = (reference_name =~ /(.*)(\.gff){0,1}(\.bed){0,1}(\.vcf){0,1}(\.gtf){0,1}/)[0][1]
      basename_fasta = (reference_name =~ /(.*)\.fasta/)[0][1]
      tagname = basename
      gzip_arg = ""
      if(gzip == params.pigz){gzip_arg = "-p ${task.cpus}"}
      cmd_gzip = "${gzip} ${gzip_arg} -c -d"
      """
      ${process_header}
      ${pigz_module}
      ${bedtools_module}
      ${cmd_gzip} ${reference_name} > ${basename_fasta}.fasta
      ${params.bedtools} getfasta -fi ${basename_fasta}.fasta -bed ${annotation_name} -fo ${basename_fasta}_split.fasta
      ${file_handle_module}
      ${file_handle_path} -r -f *_split.fasta
      """
  }
}else{
  dated_reference_names_split.into{ indexing_input }
}

process indexing {
  tag "${tagname}"
  publishDir "${index_res_path}", mode: 'copy'
  input:
    file index_name from indexing_input
    file annotation_name from dated_annotation_names_indexing
  output:
    file "*.index*" into indexing_output
    file "*_report.txt*" into indexing_log
  script:
    basename = (index_name =~ /(.*)\.fasta(\.gz){0,1}/)[0][1]
    tagname = basename
    cmd_date = "${file_handle_path} -r -f *"
    switch(mapper) {
      case "kallisto":
        """
        ${process_header}
        ${kallisto_module}
        ${params.kallisto} index -k 31 --make-unique -i ${basename}.index ${index_name} > ${basename}_kallisto_indexing_report.txt
        ${file_handle_module}
        ${cmd_date}
        """
      break
      case "bowtie2":
        gzip_arg = ""
        if(gzip == params.pigz){gzip_arg = "-p ${task.cpus}"}
        cmd_gzip = "${gzip} ${gzip_arg} -c -d"
        """
        ${process_header}
        ${pigz_module}
        ${bowtie2_module}
        ${cmd_gzip} ${index_name} > ${basename}.fasta
        ${params.bowtie2}-build --threads ${task.cpus} ${basename}.fasta ${basename}.index &> ${basename}_bowtie2_indexing_report.txt
        if grep -q "Error" ${basename}_bowtie2_indexing_report.txt; then
          exit 1
        fi
        ${file_handle_module}
        ${cmd_date}
        """
      break
      case "bowtie2+rsem":
        cmd_annotation = "--gff"
        if(annotation_name =~ /.*\.gtf/){
          cmd_annotation = "--gtf"
        }
        gzip_arg = ""
        if(gzip == params.pigz){gzip_arg = "-p ${task.cpus}"}
        cmd_gzip = "${gzip} ${gzip_arg} -c -d"
        """
        ${process_header}
        ${pigz_module}
        ${rsem_module}
        ${cmd_gzip} ${index_name} > ${basename}.fasta
        ${params.rsem}-prepare-reference -p ${task.cpus} ${rsem_parameters} ${cmd_annotation} ${annotation_name} ${basename}.fasta ${basename}.index &> ${basename}_bowtie2_rsem_indexing_report.txt
        if grep -q "Error" ${basename}_bowtie2_rsem_indexing_report.txt; then
          exit 1
        fi
        ${file_handle_module}
        ${cmd_date}
        """
      break
      default:
        """
        ${process_header}
        ${salmon_module}
        ${params.salmon} index -p ${task.cpus} -t ${index_name} -i ${basename}.index --type quasi -k 31 &> ${basename}_salmon_indexing_report.txt
        if grep -q "Error" ${basename}_salmon_indexing_report.txt; then
          exit 1
        fi
        ${file_handle_module}
        ${cmd_date}
        """
      break
    }
}

if(mapper in ["salmon", "kallisto", "bowtie2+rsem"]){
  process mapping_quantification {
    tag "${tagname}"
    publishDir "${counts_res_path}", mode: 'copy', pattern: "*.{counts,json,h5,results,cnt,model,theta,_report.txt}"
    publishDir "${bams_res_path}", mode: 'copy', pattern: "*{.bam,_report.txt}"

    input:
      file index_name from indexing_output
      file fastq_name from dated_fastq_names
    output:
      file "*.{counts,json,h5,results,cnt,model,theta,bam}" into counts_output
      file "*_report.txt" into mapping_log
    script:
      cmd_date = "${file_handle_path} -r -f"
      if (params.paired) {
        name = (fastq_name[0] =~ /(.*)\_(R){0,1}[12]\.fastq(\.gz){0,1}/)[0][1]
        tagname = name
        basename_1 = (fastq_name[0] =~ /(.*)\.fastq(\.gz){0,1}/)[0][1]
        basename_2 = (fastq_name[1] =~ /(.*)\.fastq(\.gz){0,1}/)[0][1]
        basename_index = (index_name[0] =~ /(.*\.index).*/)[0][1]
        switch(mapper) {
          case "kallisto":
            """
            ${process_header}
            ${kallisto_module}
            ${params.kallisto} quant -i ${basename_index} -t ${task.cpus} ${params.kallisto_parameters} -o ./ ${fastq_name[0]} ${fastq_name[1]} &> ${name}_kallisto_report.txt
            mv abundance.tsv ${name}.counts
            mv run_info.json ${name}_info.json
            mv abundance.h5 ${name}.h5
            if grep -q "Error" ${name}_kallisto_report.txt; then
              exit 1
            fi
            ${file_handle_module}
            ${cmd_date} *_report.txt *.counts *.json *.h5
            """
          break
          case "bowtie2+rsem":
            """
            ${process_header}
            ${rsem_module}
            ${params.rsem}-calculate-expression -p ${task.cpus} ${rsem_parameters} --paired-end ${fastq_name[0]} ${fastq_name[1]} ${basename_index} ${tagname} 2> ${name}_bowtie2_rsem_report.txt
            mv ${tagname}.stat/* .
            if grep -q "Error" ${name}_bowtie2_rsem_report.txt; then
              exit 1
            fi
            ${file_handle_module}
            ${cmd_date} "*.{results,cnt,model,theta,bam}"
            """
          break
          default:
            """
            ${process_header}
            ${salmon_module}
            ${params.salmon} quant -i ${basename_index} -p ${task.cpus} ${params.salmon_parameters} -1 ${fastq_name[0]} -2 ${fastq_name[1]} -o ${name}.counts > ${name}_salmon_report.txt
            ${file_handle_module}
            ${cmd_date} *
            """
          break
        }
      } else {
        name = (fastq_name =~ /(.*)\.fastq(\.gz){0,1}/)[0][1]
        tagname = name
        switch(mapper) {
          case "kallisto":
            kallisto_parameters = params.kallisto_parameters + " -l ${params.mean} -s ${params.sd}"
            """
            ${process_header}
            ${kallisto_module}
            ${params.kallisto} quant -i ${index_name} -t ${task.cpus} --single ${kallisto_parameters} -o ./ ${fastq_name} &> ${name}_report.txt
            mv abundance.tsv ${name}.counts
            mv run_info.json ${name}_info.json
            mv abundance.h5 ${name}.h5
            if grep -q "Error" ${name}_kallisto_report.txt; then
              exit 1
            fi
            ${file_handle_module}
            ${cmd_date} *
            """
          break
          default:
            salmon_parameters = params.salmon_parameters + " --fldMean ${params.mean} --fldSD ${params.sd}"
            """
            ${process_header}
            ${salmon_module}
            ${params.salmon} quant -i ${index_name} -p ${task.cpus} ${salmon_parameters} -r ${fastq_name} -o ${name}.counts > ${name}_report.txt
            ${file_handle_module}
            ${cmd_date} *
            """
          break
        }
      }
    }
}else{
  process mapping {
    tag "${tagname}"
    publishDir "${bams_res_path}", mode: 'copy'

    input:
      file index_name from indexing_output
      file fastq_name from dated_fastq_names
    output:
      file "*.bam" into mapping_output
      file "*_report.txt" into mapping_log
    script:
    tagname = (index_name[0] =~ /(.*)\.index.*/)[0][1]
    cmd_date = "${file_handle_path} -r -f *"
    if (params.paired) {
      name = (fastq_name[0] =~ /(.*)\_(R){0,1}[12]\.fastq(\.gz){0,1}/)[0][1]
      basename_1 = (fastq_name[0] =~ /(.*)\.fastq(\.gz){0,1}/)[0][1]
      basename_2 = (fastq_name[1] =~ /(.*)\.fastq(\.gz){0,1}/)[0][1]
      switch(mapper) {
        default:
          """
          ${process_header}
          ${bowtie2_module}
          ${params.bowtie2} ${params.bowtie2_parameters} -p ${task.cpus} -x ${tagname}.index -1 ${fastq_name[0]} -2 ${fastq_name[1]} 2> ${name}_bowtie2_report.txt | samtools view -Sb - > ${name}.bam
          if grep -q "Error" ${name}_bowtie2_report.txt; then
            exit 1
          fi
          ${file_handle_module}
          ${cmd_date}
          """
        break
      }
    } else {
      name = tagname
      switch(mapper) {
        default:
          """
          ${process_header}
          ${bowtie2_module}
          ${params.bowtie2} ${params.bowtie2_parameters} -p ${task.cpus} -x ${tagname}.index -U ${fastq_name} 2> ${name}_bowtie2_report.txt | samtools view -Sb - > ${name}.bam
          ${file_handle_module}
          ${cmd_date}
          """
        break
      }
    }
  }

  process sorting {
    tag "${tagname}"
    publishDir "${bams_res_path}", mode: 'copy'
    input:
      file bams_name from mapping_output
    output:
      file "*_sorted.bam" into sorted_mapping_output
    script:
    basename = bams_name.baseName
    tagname = basename
    """
    ${process_header}
    ${samtools_module}
    ${params.samtools} sort -@ ${task.cpus} -O BAM -o ${basename}_sorted.bam ${bams_name}
    ${file_handle_module}
    ${file_handle_path} -r -f *_sorted.bam
    """
  }

  process quantification {
    tag "${tagname}"

    echo true
    publishDir "${counts_res_path}", mode: 'copy'
    input:
      file annotation_name from dated_annotation_names_quantification
      file bams_name from sorted_mapping_output
    output:
      file "*.count*" into counts_output
    script:
    basename = bams_name.baseName
    tagname = basename
    cmd_date = "${file_handle_path} -r -f *.counts"
    if (params.paired) {
      switch(params.quantifier) {
        default:
          """
          ${process_header}
          ${htseq_module}
          ${params.htseq} -r pos ${params.htseq_parameters} --format=bam ${bams_name} ${annotation_name} &> ${basename}.count
          ${file_handle_module}
          ${cmd_date}
          """
        break
      }
    } else {
      switch(params.quantifier) {
        default:
          """
          ${process_header}
          ${htseq_module}
          ${params.htseq} -r pos ${params.htseq_parameters} --format=bam ${bams_name} ${annotation_name} &> ${basename}.count
          ${file_handle_module}
          ${cmd_date}
          """
        break
      }
    }
  }
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
