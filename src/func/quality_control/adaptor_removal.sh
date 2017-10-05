${path.params.process_header}

env

${path.cmd_unsalt_file(file)}
${path.params.cutadapt_module}
${path.cmd_adaptor_removal(task.cpu, file)}
${path.params.python2_unload_module}

${path.params.file_handle_module}
${path.cmd_date('*.{zip,html}')}
