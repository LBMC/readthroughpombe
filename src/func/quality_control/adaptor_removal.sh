${path.params.process_header}
echo "adaptor removal"

${path.cmd_unsalt_file(file)}
ls -l
${path.params.cutadapt_module}
${path.cmd_adaptor_removal(task.cpus, file)}
${path.params.python2_unload_module}

${path.params.file_handle_module}
${path.cmd_date('*_cut*.fastq.gz')}
${path.cmd_date('*_report.txt')}
ls -l
