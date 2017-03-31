
# Quality control

The first part of the pipeline is about quality control of the data.
For this we need 4 steps:
- check the quality of the fastq files
- remove (if any) the adaptor of the reads
- trim the reads to remove nucleotid of poor quality
- check the quality of the results

Those steps are executed with the `pipe/quality_control.nf` script
