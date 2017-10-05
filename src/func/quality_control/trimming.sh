${path.params.process_header}

env

${path.cmd_unsalt_file(file)}
${path.trimming_module()}
${path.cmd_trimming(task.cpu, file)}
${path.params.python2_unload_module}

${path.params.file_handle_module}
${path.cmd_date('*_trim*.fastq.gz')}
