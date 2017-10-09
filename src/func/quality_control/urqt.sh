${path.params.process_header}

env

${path.cmd_unsalt_file(file)}
${path.params.urqt_module}
${path.cmd_urqt(task.cpus, file)}

${path.params.file_handle_module}
${path.cmd_date('*_trim*.fastq.gz')}
