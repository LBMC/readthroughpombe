module load nextflow/0.25.1

# pipeline v0.1.0
nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control --fastq_files "data/fastq/*fastq.gz" --paired false --quality_threshold 20 --do_adapter_removal false -w /home/laurent/data/readthroughpombe/work/ --cpu 10

# pipeline v0.2.1
nextflow src/pipe/pipeline.nf -c src/pipeline.config -profile sge --fastq "data/fastq/*fastq.gz" --todo "fastq+cutadapt+urqt+fastqc+multiqc" -w /home/laurent/data/readthroughpombe/work/

cp data/fastq/2017_04_14_PLBD15_rrp6D_R3.fastq.gz data/fastq/2017_04_15_PLBD15_rrp6D_R3.fastq.gz
cp data/fastq/2017_04_14_PLBD4_cdc15_118_R1.fastq.gz data/fastq/2017_04_15_PLBD4_cdc15_118_R1.fastq.gz
cp data/fastq/2017_04_14_PLBD8_cut14_208_cdc15_118_R2.fastq.gz data/fastq/2017_04_15_PLBD8_cut14_208_cdc15_118_R2.fastq.gz
cp data/fastq/2017_04_14_PLBD10_rrp6D_R2.fastq.gz data/fastq/2017_04_15_PLBD10_rrp6D_R2.fastq.gz

# pipeline v0.1.0
nextflow src/pipe/quality_control.nf -c src/nextflow.config -profile quality_control_docker --fastq_files "data/fastq/2017_04_15_*.fastq.gz" --paired false --quality_threshold 20 --do_adapter_removal false -w /home/laurent/data/readthroughpombe/work/ --cpu 10

# pipeline v0.2.1
bin/nextflow src/pipe/pipeline.nf -c src/pipeline.config -profile docker --fastq "data/fastq/2017_04_15_*fastq.gz" --todo "fastq+cutadapt+urqt+fastqc+multiqc" -w /home/laurent/data/readthroughpombe/work/
