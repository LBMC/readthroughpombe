${path.params.process_header}
echo "rsem"

${path.cmd_unsalt_file(reads)}
${path.cmd_unsalt_file(index)}
ls -l

${path.params.bowtie2_module}
${path.params.rsem_module}
${path.cmd_rsem_bowtie2(task.cpus, reads, index)}

if grep -q "Error" ${tagname}_rsem_bowtie2_report.txt; then
  cat ${tagname}_rsem_bowtie2_report.txt
  exit 1
fi

${path.params.file_handle_module}
${path.cmd_date('*')}
ls -l
