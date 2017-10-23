#$ -S /bin/bash
### change logs folder
#$ -o /home/cburny/logs
#$ -e /home/cburny/logs

module load nextflow/0.25.1

# pipeline v0.1.0
nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping --fastq_files 'results/quality_control/trimming/*.fastq.gz' --mapper 'bowtie2' --paired false --reference '../../data/readthroughpombe/data/2017_09_18_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fasta.gz' --annotation '../../data/readthroughpombe/data/2017_09_19_schizosaccharomyces_pombe.chr.gtf' --quantifier 'rsem' --cpu 10 -w /home/laurent/data/readthroughpombe/work/ -resume

nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping --fastq_files "results/quality_control/trimming/2017_09_26*.fastq.gz" --mapper 'bowtie2' --paired false --reference '../../data/readthroughpombe/data/2017_09_18_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fasta.gz' --quantif false --cpu 10 -w /home/laurent/data/readthroughpombe/work/

# pipeline v0.2.1
bin/nextflow src/pipeline.nf -c src/pipeline.config -profile docker --fastq "results/quality_control/trimming/2017_09_26*.fastq.gz" --fasta '../../data/readthroughpombe/data/2017_09_18_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fasta.gz' --annot '../../data/readthroughpombe/data/2017_09_19_schizosaccharomyces_pombe.chr.gtf' --todo "fastq+cutadapt+urqt+fastqc+multiqc+bowtie2+rsem" -w /home/laurent/data/readthroughpombe/work/

bin/nextflow src/pipeline.nf -c src/pipeline.config -profile docker --fastq "results/quality_control/trimming/2017_09_26*.fastq.gz" --fasta '../../data/readthroughpombe/data/2017_09_18_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fasta.gz' --todo "fastq+cutadapt+urqt+fastqc+multiqc+bowtie2" -w /home/laurent/data/readthroughpombe/work/


bin/nextflow src/pipeline.nf -c src/pipeline.config -profile docker --fastq "results/quality_control/trimming/*.fastq.gz" --index "results/mapping/index/2017_09_27_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.index*" --annot '../../data/readthroughpombe/data/2017_09_19_schizosaccharomyces_pombe.chr.gtf' --todo "bowtie2+rsem" -w /home/laurent/data/readthroughpombe/work/
