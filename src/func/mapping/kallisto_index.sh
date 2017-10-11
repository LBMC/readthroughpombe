${path.params.process_header}

${path.cmd_unsalt_file(reference)}
ls -l

${path.params.kallisto_module}
${path.params.kallisto} index -k 31 --make-unique -i ${tagname}.index ${reference} > ${tagname}_kallisto_report.txt

${path.params.file_handle_module}
${path.cmd_date('*index*')}
${path.cmd_date('*_report.txt')}
ls -l
