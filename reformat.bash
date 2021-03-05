#!/bin/bash

perc=0.7
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

