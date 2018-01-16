require(data.table)

bed <- fread("data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr.bed", sep = "\t", h = F, stringsAsFactors = F)
bed.neg <- bed[bed$V6 %in% "-", ]
bed.pos <- bed[bed$V6 %in% "+", ]
write.table(bed.pos, "data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr.pos.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(bed.neg, "data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr.neg.bed", sep = "\t", col.names = F, row.names = F, quote = F)

# Convert in 1-based systeme coordinates
bed$V2 <- bed$V2+1
write.table(bed, "data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr1.based.bed", sep = "\t", col.names = F, row.names = F, quote = F)

# Extract + and - oriented features
bed.neg <- bed[bed$V6 %in% "-", ]
bed.pos <- bed[bed$V6 %in% "+", ]
write.table(bed.pos, "data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr1.based.pos.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(bed.neg, "data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr1.based.neg.bed", sep = "\t", col.names = F, row.names = F, quote = F)

