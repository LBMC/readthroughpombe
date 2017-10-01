${path.params.process_header}

find . -type f -name "d*" | \
sed 's/^.\\/d//g' | \
awk '{system("mv d"\$0" "\$0)}'

${path.params.cutadapt_module}
${path.cmd_adaptor_removal(task.cpu, file)}
${path.params.python2_unload_module}

${path.params.file_handle_module}
${path.cmd_date('*.{zip,html}')}
