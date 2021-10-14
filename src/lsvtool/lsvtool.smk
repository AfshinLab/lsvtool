## The pipe will start with bedpe inputs

import pkg_resources

ids, = glob_wildcards("{id,[^/]+}.bedpe") #("{id}.bedpe")
#ds = [id + "calls." for id in ids]

perc=0.8
svtype="DEL",
dist=100


xx,  = glob_wildcards("to_igv_plot/{files}.bedpe")
rule Final:
    input: 
        #expand("{filename}_filtered_sorted_merged.vcf",filename=ids),
        expand("merging/{perc}_{dist}_SVs_merged_list.bedpe",perc=perc,dist=dist),
        "to_igv_plot/common.bedpe",
        expand("to_igv_plot/{int_bedpe}.batch",int_bedpe=xx)
       

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

if len(ids)==2:
    rule intersect_bedpe:
        input: expand("merging/{perc}_{dist}_SVs_merged_list.bedpe", perc=perc,dist=dist)
        output: "to_igv_plot/common.bedpe"
        run:
            print("len IDs is :", len(ids), ids)
            shell("less {input} \| cut -f 1,2,5,12 \| grep \'1 1\' \| sed \'s\/1 1 \/both\/g\' > to_igv_plot/common.bedpe")
            shell("less {input} \| cut -f 1,2,5,12 \| grep \'1 0\' \| sed \'s\/1 0 \/first\/g\' > to_igv_plot/unique_1st.bedpe")
            shell("less {input} \| cut -f 1,2,5,12 \| grep \'0 1\' \| sed \'s\/0 1 \/second\/g\' > to_igv_plot/unique_2nd.bedpe")


    rule output_igv_batches:
        input: "to_igv_plot/{int_bedpe}.bedpe"
        output: "to_igv_plot/{int_bedpe}.batch"
        run:
            shell("bedtools igv -i {input} -slop 1000 -clps -name \| sed \'s\/collapse\/squish\/g\' > {output}")

        # bedtools igv -i $out_dir/unique_1st.bedpe -slop 1000 -clps -name | sed 's/collapse/squish/g' > $out_dir/unique_1st.batch
        # bedtools igv -i $out_dir/unique_2nd.bedpe -slop 1000 -clps -name | sed 's/collapse/squish/g' > $out_dir/unique_2nd.batch
        # bedtools igv -i $out_dir/common_both.bedpe -slop 1000 -clps -name | sed 's/collapse/squish/g' > $out_dir/common_both.batch
