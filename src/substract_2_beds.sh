##### Substract gene regions from bed peaks from MUSIC: be careful because MUSIC output peak in + orientation and our libraries are strand specific reverse
# Use bed in 0 based coordinate to do the substraction

ARRAY_fwd=(results/readthrough_analysis/Output_Music/output_cut14_wt_forward/ results/readthrough_analysis/Output_Music/output_rrp6D_wt_forward/)

# Loop over directories 
for tab in ${ARRAY_fwd[@]}; do
    echo $tab
    bin/bedtools2/bin/bedtools subtract -S -a "$tab""ERs_all.bed" -b data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr.neg.bed > "$tab""ERs_all_intersect_neg.bed" 
done 

string_to_replace_with='reverse'
for tab in ${ARRAY_fwd[@]}; do
    repl="${tab/forward/$string_to_replace_with}" &&\
    bin/bedtools2/bin/bedtools subtract -S -a "$repl""ERs_all.bed" -b data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr.pos.bed > "$repl""ERs_all_intersect_pos.bed"
done 
