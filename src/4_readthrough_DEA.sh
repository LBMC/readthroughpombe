# script to run DEA on the readthrough resutls

Rscript src/DEA.R

find results/readthrough/DEA/RT/ -name "*.csv" | while read fname;
do
  echo "$fname"
  awk -F "," 'NR >= 2 {if($7 <= 0.05) {print $1}}' $fname | sed 's/"//g' > $fname.DE.csv
  grep -F -f $fname.DE.csv results/readthrough/transcript/2018_02_14_RT_annotation_forward.bed > $fname.bed
  grep -F -f $fname.DE.csv results/readthrough/transcript/2018_02_14_RT_annotation_reverse.bed >> $fname.bed
done

find results/readthrough/DEA/T/ -name "*.csv" | while read fname;
do
  echo "$fname"
  awk -F "," 'NR >= 2 {if($7 <= 0.05) {print $1}}' $fname | sed 's/"//g' > $fname.DE.csv
  grep -F -f $fname.DE.csv results/readthrough/transcript/2018_02_14_transcript_annotation_forward.bed > $fname.bed
  grep -F -f $fname.DE.csv results/readthrough/transcript/2018_02_14_transcript_annotation_reverse.bed >> $fname.bed
done
