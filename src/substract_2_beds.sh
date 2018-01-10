##### Substract gene regions from bed peaks from MUSIC: be careful because MUSIC output peak in + orientation and our libraries are strand specific reverse
# Use bed in 0 based coordinate to do the substraction

ARRAY_fwd=(results/readthrough_analysis/Output_Music/output_cut14_wt_forward/ results/readthrough_analysis/Output_Music/output_rrp6D_wt_forward/)
ARRAY_rev=(results/readthrough_analysis/Output_Music/output_cut14_wt_reverse/ results/readthrough_analysis/Output_Music/output_rrp6D_wt_reverse/)

# Loop over directories 
for tab in ${ARRAY_fwd[@]}; do
    echo $tab
    bin/bedtools2/bin/bedtools subtract -S -a "$tab""ERs_all.bed" -b data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr.neg.bed > "$tab""ERs_all_intersect_neg.bed" 
done 

for tab in ${ARRAY_rev[@]}; do
    echo $tab
    bin/bedtools2/bin/bedtools subtract -a "$tab""ERs_all.bed" -b data/ReferenceGenomes/2017_10_25_sorted_schizosaccharomyces_pombe.chr.pos.bed > "$tab""ERs_all_intersect_pos.bed"
done 
