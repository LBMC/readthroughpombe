${path.params.process_header}

env

${path.params.pigz_module}
${path.cmd_gz(task.cpus, file)}

${path.params.file_handle_module}
${path.cmd_date("*.gz")}

find . -type f -name "*.gz" | \
sed 's/^.\\///g' | \
awk '{system("mv "\$0" d"\$0)}'
