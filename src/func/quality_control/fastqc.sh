${path.params.process_header}

env

${path.cmd_unsalt_file(file)}
${path.params.fastqc_module}
${path.cmd_fastqc(task.cpus, file)}

${path.params.file_handle_module}
${path.cmd_date('*.{zip,html}')}
