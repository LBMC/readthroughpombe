${path.params.process_header}

env

${path.cmd_unsalt_file(file)}
${path.params.cutadapt_module}
${path.cmd_adaptor_removal(task.cpus, file)}
${path.params.python2_unload_module}

ls -l
${path.params.file_handle_module}
${path.cmd_date('*_cut*.fastq.gz')}
${path.cmd_date('*_report.txt')}
