${path.params.process_header}

${path.cmd_unsalt_file(reference)}
ls -l

${path.params.pigz_module}
${path.cmd_ungz(task.cpus, reference)} > ${path.ungz_file_name(reference)}
${path.cmd_ungz(task.cpus, annotation)} > ${path.ungz_file_name(annotation)}

${path.params.bowtie2_module}
${path.params.rsem_module}
${path.cmd_rsem_bowtie2_index(task.cpus, path.ungz_file_name(reference), path.ungz_file_name(annotation))}

if grep -q "Error" ${tagname}_rsem_bowtie2_report.txt; then
  cat ${tagname}_rsem_bowtie2_report.txt
  exit 1
fi
${path.params.python2_unload_module}

${path.params.file_handle_module}
${path.cmd_date('*')}

ls -l
