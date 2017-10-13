${path.params.process_header}
echo "split fasta"

${path.cmd_unsalt_file(reference)}
${path.cmd_unsalt_file(annotation)}
ls -l

${path.params.pigz_module}
${path.cmd_ungz(task.cpus, reference)} > tmp_${path.ungz_file_name(reference)}
${path.cmd_ungz(task.cpus, annotation)} > tmp_${path.ungz_file_name(annotation)}

${path.params.bedtools_module}
${path.params.bedtools} getfasta -fullHeader -fi tmp_${path.ungz_file_name(reference)} -bed tmp_${path.ungz_file_name(annotation)} -fo ${tagname}_split.fasta

${path.cmd_gz(task.cpus, tagname + '_split.fasta')}

${path.params.file_handle_module}
${path.cmd_date('*_split*.fasta.gz')}
ls -l
