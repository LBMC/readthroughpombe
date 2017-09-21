#$ -S /bin/bash
### change logs folder
#$ -o /home/cburny/logs
#$ -e /home/cburny/logs

module load nextflow/0.25.1


nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping --fastq_files 'results/quality_control/trimming/*.fastq.gz' --mapper 'bowtie2' --paired false --reference '../../data/readthroughpombe/data/2017_09_18_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fasta.gz' --annotation '../../data/readthroughpombe/data/2017_09_19_schizosaccharomyces_pombe.chr.gtf' --quantifier 'rsem' --cpu 10 -w /home/laurent/data/readthroughpombe/work/

nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping --fastq_files 'results/quality_control/trimming/*.fastq.gz' --mapper 'bowtie2' --paired false --reference '../../data/2017_09_18_Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fasta.gz' --quantif false --cpu 10 -w /home/laurent/data/readthroughpombe/work/
