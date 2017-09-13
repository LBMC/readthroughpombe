#$ -S /bin/bash
### change logs folder
#$ -o /home/cburny/logs
#$ -e /home/cburny/logs

module load nextflow/0.25.1

nextflow src/pipe/mapping.nf -c src/nextflow.config -profile mapping_sge --fastq_files '/scratch/cburny/Input_RNASeq/*.fastq.gz' --mapper 'bowtie2' --paired false --reference '/scratch/cburny/Input_RNASeq/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fasta' --annotation '/scratch/cburny/Input_RNASeq/schizosaccharomyces_pombe.chr.gff' --quantifier 'rsem' --cpu 16
