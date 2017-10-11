${path.params.process_header}

ls -l

${path.params.file_handle_module}
${path.cmd_date("*")}
ls -l

find . -name "*" | \
grep -E "[0-9]{4}_[0-9]{2}_[0-9]{2}_" | \
sed 's/^.*\\///g' | \
awk '{system("mv "\$0" d"\$0)}'
ls -l
