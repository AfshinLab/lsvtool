#!/bin/bash

perc=0.8
tt="DEL" #"DEL"

for file in $(ls *calls.bedpe)
do
echo  ${file::${#file}-6}

python filter_lsv.py -f $file -t $tt -q $perc -d 100 -bl ../hg38_black_list.bed,../hg38_gap.bed,centromers.hg38.bed #DUP INV

done


ls *filtered_sorted_merged.vcf > files_to_merge.txt

SURVIVOR merge files_to_merge.txt 5000 1 0 0 0 10000 ${perc}_${tt}_naibrs_merged.vcf

Rscript plot.R -i ${perc}_${tt}_naibrs_merged.vcf -t $tt

rm *log


#if two files only use:
less ${perc}_${tt}_naibrs_merged_list.bedpe | cut -f 1,2,5,12 | grep "1 0\|1 1" | sed 's/1 0 /first/g'  | sed 's/1 1 /both/g' > All_from_1st.bedpe
less ${perc}_${tt}_naibrs_merged_list.bedpe | cut -f 1,2,5,12 | grep "0 1\|1 1" | sed 's/0 1 /second/g' | sed 's/1 1 /both/g' > All_from_2nd.bedpe

bedtools igv -i All_from_1st.bedpe -slop 1000 -clps -name | sed 's/collapse/squish/g' > all_1st.batch
bedtools igv -i All_from_2nd.bedpe -slop 1000 -clps -name | sed 's/collapse/squish/g' > all_2nd.batch