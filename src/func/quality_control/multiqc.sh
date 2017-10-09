${path.params.process_header}

env

${path.params.multiqc_module}
${path.params.multiqc} -f .
${path.params.python2_unload_module}

${path.params.file_handle_module}
ls -l
${path.cmd_date('multiqc_*')}
ls -l
