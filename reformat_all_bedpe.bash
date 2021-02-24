#!/bin/bash


for file in $(ls *calls.bedpe)
do
echo  ${file::${#file}-6}

python filter_lsv.py -f $file -t DEL -q 0.50 -d 10 -bl ../hg38_black_list.bed,../hg38_gap.bed #DUP INV

done


ls *sorted_merged.vcf > files_to_merge.txt

SURVIVOR merge files_to_merge.txt 5000 1 0 0 0 10000 4_merged.vcf
Rscript plot.R -i 4_merged.vcf

rm *log
