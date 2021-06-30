#!/bin/bash

perc=0.8
tt="DEL" #"DEL"

out_dir="output"$perc$tt

for file in $(ls *calls.bedpe)
do
echo  ${file::${#file}-6}

python filter_lsv.py -f $file -t $tt -q $perc -d 100 -o $out_dir #-bl ../hg38_black_list.bed,../hg38_gap.bed,centromers.hg38.bed #DUP INV

done


ls $out_dir/*filtered_sorted_merged.vcf > $out_dir/files_to_merge.txt

SURVIVOR merge $out_dir/files_to_merge.txt 5000 1 0 0 0 10000 $out_dir/${perc}_${tt}_naibrs_merged.vcf

Rscript plot.R -i $out_dir/${perc}_${tt}_naibrs_merged.vcf -t $tt 

rm $out_dir/*log


#if two files only use:
less $out_dir/${perc}_${tt}_naibrs_merged_list.bedpe | cut -f 1,2,5,12 | grep "1 0\|1 1" | sed 's/1 0 /first/g'  | sed 's/1 1 /both/g' > $out_dir/All_from_1st.bedpe
less $out_dir/${perc}_${tt}_naibrs_merged_list.bedpe | cut -f 1,2,5,12 | grep "0 1\|1 1" | sed 's/0 1 /second/g' | sed 's/1 1 /both/g' > $out_dir/All_from_2nd.bedpe

bedtools igv -i $out_dir/All_from_1st.bedpe -slop 1000 -clps -name | sed 's/collapse/squish/g' > $out_dir/All_1st.batch
bedtools igv -i $out_dir/All_from_2nd.bedpe -slop 1000 -clps -name | sed 's/collapse/squish/g' > $out_dir/All_2nd.batch
