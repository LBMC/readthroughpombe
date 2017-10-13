${path.params.process_header}

ls -l
${path.cmd_unsalt_file(index)}
${path.cmd_unsalt_file(reads)}
ls -l

${path.params.kallisto_module}
${path.cmd_kallisto(task.cpus, index, reads)}

mv abundance.tsv ${tagname}.tsv
mv run_info.json ${tagname}_info.json
mv abundance.h5 ${tagname}.h5
if grep -q "Error" ${tagname}_kallisto_report.txt; then
  exit 1
fi

${path.params.file_handle_module}
${path.cmd_date('*')}
ls -l
