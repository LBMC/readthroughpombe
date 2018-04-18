# script to run DEA on the readthrough resutls

Rscript src/DEA.R
Rscript src/DEG_set.R

find results/readthrough/DEA/RT/ -name "*.csv" | while read fname;
do
  echo "$fname"
  awk -F "," 'NR >= 2 {if($7 <= 0.05) {print $1}}' $fname | sed 's/"//g' > $fname.DE.csv
  grep -F -f $fname.DE.csv results/readthrough/transcript/2018_02_14_RT_annotation_forward.bed > $fname.bed
  grep -F -f $fname.DE.csv results/readthrough/transcript/2018_02_14_RT_annotation_reverse.bed >> $fname.bed
  awk -F "," 'NR >= 2 {print $1}' $fname | sed 's/"//g' > $fname.all.csv
  grep -F -f $fname.all.csv results/readthrough/transcript/2018_02_14_RT_annotation_forward.bed > $fname.all.bed
  grep -F -f $fname.all.csv results/readthrough/transcript/2018_02_14_RT_annotation_reverse.bed >> $fname.all.bed
done

find results/readthrough/DEA/T/ -name "*.csv" | while read fname;
do
  echo "$fname"
  awk -F "," 'NR >= 2 {if($7 <= 0.05) {print $1}}' $fname | sed 's/"//g' > $fname.DE.csv
  grep -F -f $fname.DE.csv results/readthrough/transcript/2018_02_14_transcript_annotation_forward.bed > $fname.bed
  grep -F -f $fname.DE.csv results/readthrough/transcript/2018_02_14_transcript_annotation_reverse.bed >> $fname.bed
  awk -F "," 'NR >= 2 {print $1}' $fname | sed 's/"//g' > $fname.all.csv
  grep -F -f $fname.all.csv results/readthrough/transcript/2018_02_14_transcript_annotation_forward.bed > $fname.all.bed
  grep -F -f $fname.all.csv results/readthrough/transcript/2018_02_14_transcript_annotation_reverse.bed >> $fname.all.bed
done
