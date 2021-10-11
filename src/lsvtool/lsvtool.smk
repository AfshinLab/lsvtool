## The pipe will start with bedpe inputs

import pkg_resources

ids, = glob_wildcards("{id}.bedpe")

perc=0.8
svtype="DEL",
dist=100

scattergather:
    split=8

rule Final:
    input: 
        #expand("{filename}_filtered_sorted_merged.vcf",filename=ids),
        expand("{perc}_{dist}_SVs_merged_overlap.png",perc=perc,dist=dist),
        #"common_both.bedpe"
        #"common_both.bedpe"

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

        print("\n\n\n",len(ids),ids,"\n\n\n")
        shell("lsvtool filter_lsv  -f {input} -t {svtype} -q {perc} -d {dist} -bl {defaultBlacklists} -o filtered_inputs")

rule merge_lsv_files:
    input: expand("filtered_inputs/{filename}_filtered_sorted_merged.vcf", filename=ids)
    output: 
        "{perc}_{dist}_SVs_merged.vcf"
    params: 
        merge = "files_to_merge.txt"
    shell:
        "cd filtered_inputs && ls *_filtered_sorted_merged.vcf > {params.merge}"
        " && SURVIVOR merge {params.merge} 5000 1 0 0 0 10000 ../{output}"


rule plot_intersection:
    input: "{perc}_{dist}_SVs_merged.vcf"
    output:"{perc}_{dist}_SVs_merged_overlap.png"
    params: pkg_resources.resource_filename("lsvtool", "cli/plot.R")
    shell: 
        "Rscript {params} -i {input}"


# rule output_igv_batches:
#     input: expand("{perc}_{dist}_SVs_merged.vcf",perc=perc,dist=dist)
#     output: "common_both.bedpe"
#     run:
#         twofiles = ["less {input} \| cut -f 1,2,5,12 \| grep \'1 0\' \| sed \'s\/1 0 \/first\/g\' > common_both.bedpe"]
#         #if two files only use:
#         print("len IDs is :", len(ids), IDs)
#         if len(ids)==2:
#             shell(commands[0]) # && \\
                     #less {input} | cut -f 1,2,5,12 | grep '0 1' | sed 's/0 1 /second/g' > unique_2nd.bedpe && \\ 
                     #less {input} | cut -f 1,2,5,12 | grep '1 1' | sed 's/1 1 /both/g' > common_both.bedpe """)


# less $out_dir/${perc}_${tt}_SVs_merged_list.bedpe | cut -f 1,2,5,12 | grep "1 0" | sed 's/1 0 /first/g'   > $out_dir/unique_1st.bedpe
# less $out_dir/${perc}_${tt}_SVs_merged_list.bedpe | cut -f 1,2,5,12 | grep "0 1" | sed 's/0 1 /second/g'  > $out_dir/unique_2nd.bedpe
# less $out_dir/${perc}_${tt}_SVs_merged_list.bedpe | cut -f 1,2,5,12 | grep "1 1" | sed 's/1 1 /both/g'    > $out_dir/common_both.bedpe

# bedtools igv -i $out_dir/unique_1st.bedpe -slop 1000 -clps -name | sed 's/collapse/squish/g' > $out_dir/unique_1st.batch
# bedtools igv -i $out_dir/unique_2nd.bedpe -slop 1000 -clps -name | sed 's/collapse/squish/g' > $out_dir/unique_2nd.batch
# bedtools igv -i $out_dir/common_both.bedpe -slop 1000 -clps -name | sed 's/collapse/squish/g' > $out_dir/common_both.batch
