${path.params.process_header}

env

ls -l
${path.params.pigz_module}
${path.cmd_gz(task.cpus, file)}
ls -l

${path.params.file_handle_module}
${path.cmd_date("*.gz")}
ls -l

find . -name "*.gz" | \
sed 's/^.*\\///g' | \
awk '{system("mv "\$0" d"\$0)}'
ls -l
