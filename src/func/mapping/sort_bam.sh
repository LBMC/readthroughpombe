${path.params.process_header}
echo "sort bam"

${path.cmd_unsalt_file(bam)}
ls -l

${path.params.samtools_module}
${path.params.samtools} sort -@ ${task.cpus} -O BAM -o ${tagname}_sorted.bam ${path.unsalt_file_name(bam)}

${path.params.file_handle_module}
${path.cmd_date('*_sorted.bam')}
ls -l
