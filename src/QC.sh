module load nextflow/0.25.1

nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_sge --fastq_files "data/fastq/*fastq.gz" --paired false --quality_threshold 20 --do_adapter_removal false

