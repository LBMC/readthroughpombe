${path.params.process_header}

find . -type f -name "d*" | \
sed 's/^.\\/d//g' | \
awk '{system("mv d"\$0" "\$0)}'

${path.trimming_module()}
${path.cmd_trimming(task.cpu, file)}
${path.params.python2_unload_module}

${path.params.file_handle_module}
${path.cmd_date('*_trim*.fastq.gz')}
