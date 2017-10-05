${process_header}
${kallisto_module}
${params.kallisto} quant -i ${basename_index} -t ${task.cpus} ${params.kallisto_parameters} -o ./ ${fastq_name[0]} ${fastq_name[1]} &> ${name}_kallisto_report.txt
mv abundance.tsv ${name}.counts
mv run_info.json ${name}_info.json
mv abundance.h5 ${name}.h5
if grep -q "Error" ${name}_kallisto_report.txt; then
  exit 1
fi
${file_handle_module}
${cmd_date} *_report.txt *.counts *.json *.h5
