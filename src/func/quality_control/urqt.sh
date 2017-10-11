${path.params.process_header}
echo "UrQt"

${path.cmd_unsalt_file(file)}
ls -l

${path.params.urqt_module}
${path.cmd_urqt(task.cpus, file)}

${path.params.file_handle_module}
${path.cmd_date('*_trim*.fastq.gz')}
ls -l
