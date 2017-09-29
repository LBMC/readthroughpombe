${path.params.process_header}

${path.params.pigz_module}
${path.cmd_gz(task.cpus, file_S)}

${path.params.file_handle_module}
${path.cmd_date(file_S + ".gz")}

find . -type f -name "*.fastq.gz" | \
sed 's/^.\\///g' | \
awk '{system("mv "\$0" d"\$0)}'
