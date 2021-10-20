## The pipe will start with bedpe inputs

import pkg_resources
import glob
import numpy as np

ids, = glob_wildcards("{id,[^/]+}.bedpe") 

perc=0.8
svtype="DEL",
dist=1000

rule Final:
    input: 
        #expand("{filename}_filtered_sorted_merged.vcf",filename=ids),
        expand("merging/{perc}_{dist}_SVs_merged_list.bedpe",perc=perc,dist=dist),
        "to_igv_plot/common.batch"
       

rule filter_lsv_files:
    input:
        "{filename}.bedpe"
    output: 
        "filtered_inputs/{filename}_filtered_sorted_merged.vcf"
    run:
        gap_bed = pkg_resources.resource_filename("lsvtool", "refs/hg38_gap.bed")
        centromers_bed = pkg_resources.resource_filename("lsvtool", "refs/hg38_centromers.bed")
        blacklist_bed = pkg_resources.resource_filename("lsvtool", "refs/hg38_black_list.bed")
        defaultBlacklists = f"{blacklist_bed},{gap_bed},{centromers_bed}"
        #print("\n\n\n",len(ids),ids,"\n\n\n")
        shell("lsvtool filter_lsv  -f {input} -t {svtype} -q {perc} -d {dist} -bl {defaultBlacklists} -o filtered_inputs")

rule merge_lsv_files:
    input: expand("filtered_inputs/{filename}_filtered_sorted_merged.vcf", filename=ids)
    output: 
        "merging/{perc}_{dist}_SVs_merged.vcf"
    params: 
        merge = "files_to_merge.txt"
    shell:
        "cd filtered_inputs && ls *_filtered_sorted_merged.vcf > {params.merge}"
        " && SURVIVOR merge {params.merge} 5000 1 0 0 0 10000 ../{output}"


rule plot_intersection:
    input: "merging/{perc}_{dist}_SVs_merged.vcf"
    output:"merging/{perc}_{dist}_SVs_merged_list.bedpe"
    params: pkg_resources.resource_filename("lsvtool", "cli/plot.R")
    shell: 
        "Rscript {params} -i {input}"


rule intersect_bedpe:
    input: expand("merging/{perc}_{dist}_SVs_merged_list.bedpe", perc=perc,dist=dist)
    output: "to_igv_plot/{files}.bedpe"
    run:
        n = len(ids)
        xx = np.zeros((n,n))
        for i in range (0,n):
            xx[i,i]=1
        xx = np.vstack( [xx,np.ones((n))] )
        combinations = [' '.join(c for c in str(xx[line]) if c.isdigit()) for line in range (0,n+1)]
        
        for i in combinations :
            lable=i.replace(" ", "")
            if not "0" in lable :
                lable="common"
            try:
                shell(f"less {input} | cut -f 1,2,5,12 | grep \'{i}\' | sed \'s/{i} /{lable}/g\' > to_igv_plot/{lable}.bedpe")
            except:
                continue

rule output_igv_batches:
    input: "to_igv_plot/{file}.bedpe"
    output:  "to_igv_plot/{file}.batch"
    run:
        for file in glob.glob("to_igv_plot/*bedpe"):
            fileout = file.rsplit('.', 1)[0] + '.batch' 
            shell("bedtools igv -i {file} -slop 1000 -clps -name | sed \'s/collapse/squish/g\' > {fileout}")
