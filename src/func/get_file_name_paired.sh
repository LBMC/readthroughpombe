${path.params.process_header}

${path.params.pigz_module}
${path.cmd_gz(task.cpus, file_R1)}
${path.cmd_gz(task.cpus, file_R2)}

${path.params.file_handle_module}
${path.cmd_date(file_R1 + ".gz")}
${path.cmd_date(file_R2 + ".gz")}

find . -type f -name "*.fastq.gz" | \
sed 's/^.\\///g' | \
awk '{system("mv "\$0" d"\$0)}'
