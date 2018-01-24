# script to build hsp9 example to work on

mkdir -p data/example/

# hsp9 position in the file data/2017_09_19_schizosaccharomyces_pombe.chr.gff3
# I:5319957-5320676
# we take the region with 100kb on both side of hsp9

# get annotation
grep -e "I" data/2017_09_19_schizosaccharomyces_pombe.chr.gff3 | \
awk -v FS='\t' -v OFS='\t' '{if($4 >= 5219957 && $5 <= 5420676){$4=$4-5219957; $5=$5-5219957; print $0}}' > data/example/example.gff3

grep -e "I" data/2017_09_19_schizosaccharomyces_pombe.chr.gtf | \
awk -v FS='\t' -v OFS='\t' '{if($4 >= 5219957 && $5 <= 5420676){$4=$4-5219957; $5=$5-5219957; print $0}}' > data/example/example.gtf

# get sequence
echo "I\t5219957\t5420676\ttoextract" > data/example/getfasta.bed
bedtools getfasta -fi data/Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa -bed data/example/getfasta.bed -fo data/example/example.fasta
sed -i 's/>I:5219957-5420676/>I/g' data/example/example.fasta

#get mapping
find results/mapping/mapping/ -name "2018_01_19_*_rev_sort.bam" | \
sed 's/results\/mapping\/mapping\/2018_01_19_//g' | \
awk '{system("samtools index results/mapping/mapping/2018_01_19_"$1"; samtools view -@ 11 -b results/mapping/mapping/2018_01_19_"$1" -o data/example/example_"$1" I:5219958-5420675")}'

find data/example/*.bam | \
sed 's/bam//g' | \
awk '{system("bedtools bamtofastq -i "$1"bam -fq "$1"fastq")}'

#we have everything in the data/example folder

################################################################################
################################## Mapping #####################################

bin/nextflow src/pipeline.nf -c src/pipeline.config -profile docker --fastq "data/example/*.fastq" --fasta 'data/example/example.fasta' --todo "bowtie2"

################################################################################
########################### readthrough detection ##############################

bin/nextflow src/readthrough.nf -c src/pipe/conf/readthrough_docker.config --name "rrp6D" --bams "data/example/bams_rrp6D.csv" --annotation "data/example/example.gff3" --index "results/mapping/indexing/2018_01_18_example.index*" --genome "data/example/example.fasta" --reads_size 50 --frag_size 363

bin/nextflow src/readthrough.nf -c src/pipe/conf/readthrough_docker.config --name "cut14" --bams "data/example/bams_cut14.csv" --annotation "data/example/example.gff3" --index "results/mapping/indexing/2018_01_18_example.index*" --genome "data/example/example.fasta" --reads_size 50 --frag_size 363 -resume

bin/nextflow src/readthrough.nf -c src/pipe/conf/readthrough_docker.config --name "cut14_cdc15" --bams "data/example/bams_cut14_cdc15.csv" --annotation "data/example/example.gff3" --index "results/mapping/indexing/2018_01_18_example.index*" --genome "data/example/example.fasta" --reads_size 50 --frag_size 363 -resume
