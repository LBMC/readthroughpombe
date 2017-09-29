${path.params.process_header}

find . -type f -name "d*" | \
sed 's/^.\\/d//g' | \
awk '{system("mv d"\$0" "\$0)}'

${path.params.fastqc_module}
${path.cmd_fastqc(task.cpu, file)}

${path.params.file_handle_module}
${path.cmd_date('*.{zip,html}')}
