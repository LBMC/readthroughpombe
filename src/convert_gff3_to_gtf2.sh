# Convert gff3 annotation to gtf2 for RSEM input 

bin/gffread -E $1 -T -o- > $2
