${path.params.process_header}

${path.cmd_unsalt_file(reference)}
ls -l

${path.params.pigz_module}
${path.cmd_ungz(task.cpus, reference)} > tmp_${path.ungz_file_name(reference)}

${path.params.bowtie2_module}
${path.params.bowtie2}-build --threads ${task.cpus} tmp_${path.ungz_file_name(reference)} ${tagname}.index &> ${tagname}_bowtie2_report.txt

if grep -q "Error" ${tagname}_bowtie2_report.txt; then
  exit 1
fi

${path.params.file_handle_module}
${path.cmd_date('*index*')}
${path.cmd_date('*_report.txt')}
ls -l
