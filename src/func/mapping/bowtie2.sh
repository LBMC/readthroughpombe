${path.params.process_header}

ls -l
${path.cmd_unsalt_file(index)}
${path.cmd_unsalt_file(reads)}
ls -l

${path.params.bowtie2_module}
${path.cmd_bowtie2(task.cpus, index, reads)}

if grep -q "Error" ${tagname}_bowtie2_report.txt; then
  exit 1
fi

${path.params.file_handle_module}
${path.cmd_date('*')}
ls -l
