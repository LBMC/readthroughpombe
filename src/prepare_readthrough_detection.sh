##### Sort annotation
sort -k1,1 -k4,4n data/ReferenceGenomes/schizosaccharomyces_pombe.chr.gff3 > data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr.gff3 &&\

##### Convert gff annotation to bed 
bash src/convert_gff3_to_bed.sh data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr.gff3 data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr.bed &&\

##### Extract + and - features in 1 based coordinates
Rscript src/format_bed.R &&\

##### Index sorted BAM
for filename in results/mapping/mapped/cut14-208/*.bam; do
    samtools index "$filename" 
done

for filename in results/mapping/mapped/wt/*.bam; do
    samtools index "$filename" 
done

for filename in results/mapping/mapped/cdc15-118/*.bam; do
    samtools index "$filename" 
done

for filename in results/mapping/mapped/rrp6D/*.bam; do
    samtools index "$filename" 
done

for filename in results/mapping/mapped/cut14-208_cdc15-118/*.bam; do
    samtools index "$filename" 
done

##### Consider sense and antisense reads
for filename in results/mapping/mapped/cut14-208/*sort.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb -F 0x10 "$filename" > "$res""_forward.bam" &&\
    samtools view -hb -f 0x10 "$filename" > "$res""_reverse.bam" 
done

for filename in results/mapping/mapped/wt/*sort.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb -F 0x10 "$filename" > "$res""_forward.bam" &&\
    samtools view -hb -f 0x10 "$filename" > "$res""_reverse.bam" 
done

for filename in results/mapping/mapped/cdc15-118/*sort.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb -F 0x10 "$filename" > "$res""_forward.bam" &&\
    samtools view -hb -f 0x10 "$filename" > "$res""_reverse.bam" 
done

for filename in results/mapping/mapped/rrp6D/*sort.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb -F 0x10 "$filename" > "$res""_forward.bam" &&\
    samtools view -hb -f 0x10 "$filename" > "$res""_reverse.bam" 
done

for filename in results/mapping/mapped/cut14-208_cdc15-118/*sort.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb -F 0x10 "$filename" > "$res""_forward.bam" &&\
    samtools view -hb -f 0x10 "$filename" > "$res""_reverse.bam" 
done

##### Filter for BAM selected regions outside annotation
for filename in results/mapping/mapped/cut14-208/*sort_forward.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb "$filename" -L data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr1.based.neg.bed -U "$res""_without.bam" > /dev/null
    samtools index "$res""_without.bam"
done 

for filename in results/mapping/mapped/wt/*sort_forward.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb "$filename" -L data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr1.based.neg.bed -U "$res""_without.bam" > /dev/null
    samtools index "$res""_without.bam"
done 

for filename in results/mapping/mapped/cdc15-118/*sort_forward.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb "$filename" -L data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr1.based.neg.bed -U "$res""_without.bam" > /dev/null
    samtools index "$res""_without.bam"
done 

for filename in results/mapping/mapped/rrp6D/*sort_forward.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb "$filename" -L data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr1.based.neg.bed -U "$res""_without.bam" > /dev/null
    samtools index "$res""_without.bam"
done 

for filename in results/mapping/mapped/cut14-208_cdc15-118/*sort_forward.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb "$filename" -L data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr1.based.neg.bed -U "$res""_without.bam" > /dev/null
    samtools index "$res""_without.bam"
done 

for filename in results/mapping/mapped/cut14-208/*sort_reverse.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb "$filename" -L data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr1.based.pos.bed -U "$res""_without.bam" > /dev/null
    samtools index "$res""_without.bam"
done 

for filename in results/mapping/mapped/wt/*sort_reverse.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb "$filename" -L data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr1.based.pos.bed -U "$res""_without.bam" > /dev/null
    samtools index "$res""_without.bam"
done 

for filename in results/mapping/mapped/cdc15-118/*sort_reverse.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb "$filename" -L data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr1.based.pos.bed -U "$res""_without.bam" > /dev/null
    samtools index "$res""_without.bam"
done 

for filename in results/mapping/mapped/rrp6D/*sort_reverse.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb "$filename" -L data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr1.based.pos.bed -U "$res""_without.bam" > /dev/null
    samtools index "$res""_without.bam"
done 

for filename in results/mapping/mapped/cut14-208_cdc15-118/*sort_reverse.bam; do
    res="$(echo $filename | cut -d'.' -f1)" &&\
    samtools view -hb "$filename" -L data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr1.based.pos.bed -U "$res""_without.bam" > /dev/null
    samtools index "$res""_without.bam"
done 
