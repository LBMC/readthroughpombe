${path.params.process_header}
echo "HTSeq"

${path.cmd_unsalt_file(bam)}
${path.cmd_unsalt_file(annotation)}
ls -l

${path.params.htseq_module}
${path.params.htseq} -r pos ${path.params.htseq_parameters} --format=bam ${path.unsalt_file_name(bam)} ${path.unsalt_file_name(annotation)} &> ${tagname}.count

${path.params.python2_unload_module}
${path.params.file_handle_module}
${path.cmd_date('*.count')}
ls -l
