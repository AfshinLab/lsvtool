import glob
from itertools import product

from lsvtool.utils import check_path

configfile: "parameters.config"


svtype = config.get("svtype", "ALL")
dist_merge = float(config.get("dist_merge", 1000))

# 20kbp is used by longranger (see https://github.com/10XGenomics/longranger/blob/e2a3143b3956af6290fd4ba08e09f76985293685/mro/_combined_sv_caller.mro#L91)
segdups_dist = int(config.get("dist_segdups", 20_000))  

minlength = int(config.get("minlength", 50))
maxlength = int(config.get("maxlength", 1_000_000_000))
BL = check_path(config.get('BL'))
segdups = check_path(config.get('segdups'))

bench_bed = check_path(config.get('bench_bed'))
bench_vcf = check_path(config.get('bench_vcf'))
igv_batches = bool(config.get('igv_batches', True))

filters = config.get("filters", None)
collapse_args = config.get("collapse_args", "")
bench_args = config.get("bench_args", "")

# Get samples from VCF names
samples = set(glob_wildcards("{id,[^/]+}.vcf.gz").id)  # Gzipped
samples |= set(glob_wildcards("{id,[^/]+}.vcf").id) # Unzipped
samples = list(sorted(samples))

fig_formats = ["png", "svg"]

final_input = [
    expand("3_filtered/{filename}.final.bedpe", filename=samples),
    expand("3_filtered/{filename}.final.vcf.gz.tbi", filename=samples),
    "4_merged/merged.vcf",
    "4_merged/merged.jl",
    expand("4_merged/merged_comp_mat_heatmap.{ext}", ext=fig_formats),
    expand("4_merged/merged_comp_mat_jaccard_heatmap.{ext}", ext=fig_formats),
    expand("4_merged/merged_venn.{ext}", ext=fig_formats),   
]


if bench_vcf:
    final_input.extend([
        expand("bench/{filename}.jl", filename=samples)
    ])

if igv_batches:
    expected_files = ["".join(p) for p in product("01", repeat=len(samples))]
    expected_files.remove("1"*len(samples))  # Called 'common'
    expected_files.remove("0"*len(samples))  # No calls
    expected_files.append("common")

    final_input.extend([
        expand("5_to_igv_plot/{files}.batch", files=expected_files),
        "5_to_igv_plot/README.md"
    ])


rule final:
    input: 
        final_input


rule bgzip:
    input:
        vcf = "{file}.vcf"
    output:
        vcf = "{file}.vcf.gz"
    shell:
        "bcftools sort {input.vcf} | bgzip -c > {output.vcf}"


rule tabix:
    input:
        vcf = "{file}.vcf.gz"
    output:
        vcf = "{file}.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input.vcf}"


ruleorder: select > bgzip


rule select:
    """Select variants of correct type that passed filters"""
    input:
        vcf = "{filename}.vcf.gz"
    output:
        vcf = "1_selected/{filename}.vcf.gz"
    params: 
        select = '' if svtype == "ALL" else f"-i 'INFO/SVTYPE == \"{svtype}\"'",
        filters = '' if filters is None else f"-f '{filters}'"
    shell:
        "bcftools view {params.filters} {params.select} {input.vcf} | bgzip -c > {output.vcf}"


rule collapse:
    """Collapse variants of similar size and breakpoints"""
    input:
        vcf = "1_selected/{filename}.vcf.gz",
        index = "1_selected/{filename}.vcf.gz.tbi",
    output:
        vcf = "2_collapsed/{filename}.vcf",
        removed = "2_collapsed/{filename}.removed.vcf"
    log: "2_collapsed/{filename}.vcf.log"
    params:
        sizemin = int(minlength * 0.8),
        sizemax = int(maxlength / 0.8)
    run:
        command = (
            "truvari collapse"
            " --input {input.vcf}"
            " --collapsed-output {output.removed}"
            " --keep *"
            " {collapse_args}"
            " --sizemin {params.sizemin}"
            " --sizemax {params.sizemax}"
            " 2> {log}"
            " |"
            " bcftools sort -"
            " > {output.vcf}"
            " 2>> {log}")
        try: # Try taking the SV with max quality if quality provided
            shell(command.replace("*", "maxqual"))
        except:
            shell(command.replace("*", "first"))


rule filter_blacklist:
    """Filter SVs against blacklist and for size"""
    input:
        vcf = "2_collapsed/{filename}.vcf"
    output:
        vcf = "3_filtered/{filename}.vcf"
    log: "3_filtered/{filename}.vcf.log"
    params:
        blacklist = BL if BL is not None else "NA"
    shell:
        "SURVIVOR filter"
        " {input.vcf}"
        " {params.blacklist}"
        " {minlength}"
        " {maxlength}"
        " -1"
        " -1"
        " {output.vcf}"
        " > {log}"


rule filter_segdups:
    """Filter SVs for overlap with segmental duplications if provided"""
    input:
        vcf = "3_filtered/{filename}.vcf"
    output:
        vcf = "3_filtered/{filename}.final.vcf"
    log: "3_filtered/{filename}.final.vcf.log"
    run:
        if segdups:
            shell(
                "lsvtool intersect_bedpe"
                " {input.vcf}"
                " {segdups}"
                " -o {output.vcf}"
                " -d {segdups_dist}" 
                " 2> {log}"
            )
        else:
            shell("ln -s $PWD/{input.vcf} $PWD/{output.vcf}")


rule vcftobedpe:
    """Convert from VCF to BEDPE"""
    # NB: Comman 'SURVIVOR vcftobed' acctually outputs BEDPE
    input:
        vcf = "{filename}.vcf"
    output:
        bedpe = "{filename}.bedpe"
    shell:
        "SURVIVOR vcftobed"
        " {input.vcf}"
        " -1"
        " -1"
        " {output.bedpe}"


rule merge_svs:
    """Merge SV between multiple VCF to find common variants"""
    input: 
        vcfs = expand("3_filtered/{filename}.final.vcf", filename=samples)
    output: 
        vcf = temp("4_merged/merged_no_names.vcf"),
        list = "4_merged/files_merged.list"
    log: "4_merged/merged_no_names.vcf.log"
    run:
        with open(output.list, "w") as f:
            f.writelines("\n".join(input.vcfs) + "\n")

        shell("SURVIVOR merge"
              " {output.list}"
              " {dist_merge}"
              " 1"
              " 0"
              " 0"
              " 0"
              " {minlength}"
              " {output.vcf}"
              " > {log}")


rule reheader:
    """Add correct sample names to files"""
    input:
        vcf = "4_merged/merged_no_names.vcf"
    output:
        vcf = "4_merged/merged.vcf",
        list = "4_merged/sample_names.list"
    run:
        with open(output.list, "w") as f:
            f.writelines("\n".join(samples) + "\n")
        shell("bcftools reheader -s {output.list} {input.vcf} > {output.vcf}")


rule gencomp:
    """Compare SVs between sample and output comparison matrix"""
    input:
        vcf = "4_merged/merged.vcf"
    output:
        txt = "4_merged/merged_comp_mat.txt"
    shell:
        "SURVIVOR genComp {input.vcf} 0 {output.txt}"


rule plot_heatmap:
    """Plot matrix as heatmap"""
    input:
        "4_merged/merged_comp_mat.txt",
        "4_merged/sample_names.list"
    output:
        "4_merged/merged_comp_mat_heatmap.{ext}",
        "4_merged/merged_comp_mat_jaccard_heatmap.{ext}",
    script:
        "scripts/plot_heatmap.R"


rule intersection:
    """For each SV list the samples it appears for"""
    # Code from https://github.com/fritzsedlazeck/SURVIVOR/wiki#plotting-the-comparison-of-multiple-input-vcf-files-after-merging
    input:
        "4_merged/merged.vcf"
    output:
        "4_merged/merged_intersection.txt"
    shell:
        "perl -ne \'print \"$1\\n\" if /SUPP_VEC=([^,;]+)/\' {input} | sed -e \'s/\\(.\\)/\\1 /g\' > {output}"


rule plot_venn:
    """Plot venn-diagram of SVs"""
    input:
        "4_merged/merged_intersection.txt",
        "4_merged/sample_names.list"
    output:
        "4_merged/merged_venn.{ext}",
    params:
        svtype = svtype
    script:
        "scripts/plot_venn.R"


ruleorder: intersect_bedpe > vcftobedpe


rule intersect_bedpe:
    # Attach intersection info to BEDPE
    input: 
        bedpe = "4_merged/merged.bedpe",
        intersection = "4_merged/merged_intersection.txt",
    output: 
        joint_bedpe = "4_merged/merged_intersection.bedpe",
    shell:
        "paste {input.bedpe} {input.intersection} > {output.joint_bedpe}"


rule get_igv_regions:
    input:
        bedpe = "4_merged/merged_intersection.bedpe",
    output:
        bed = touch("5_to_igv_plot/{files}.bed")
    run:
        label=wildcards.files
        i = " ".join(list(label))
        if label == "common":
            i = " ".join(list("1"*len(samples)))
        try:
            shell(f"cut -f 1,2,5,12 {input.bedpe} | grep \'{i} $\' | sed \'s/{i} /{label}/g\' > {output.bed}")
        except Exception:
            # Runtime error raised if `i` is missing from input
            pass


rule output_igv_batches:
    input: 
        bed = "5_to_igv_plot/{file}.bed"
    output: 
        batch = "5_to_igv_plot/{file}.batch"
    shell:
        "bedtools igv"
        " -i {input.bed}"
        " -slop 1000"
        " -clps"
        " -name"
        " -path {wildcards.file}"
        " |"
        " sed \'s/collapse/squish/g\'"
        " > {output.batch}"


rule to_igv_plot_readme:
    output:
        "5_to_igv_plot/README.md"
    run:
        with open(output[0], "w") as fout:
            print("""
This directory contains to BED and BATCH files useful to visualize and compare SV found 
in different samples.

BED files (`*.bed`) contain regions for SV calls. 

BATCH files (`*.batch`) contain scripts to generate images of said regions within a IGV session.

See table below for which files contain which sets.
""", file=fout)
            wcol = max(map(len, samples)) + 2
            wname = max(max(map(len, expected_files)) + 1, 7)
            true_str = "X".center(wcol, " ")
            false_str = " "*wcol
            print("Prefix".ljust(wname, " "), *[i.center(wcol, " ") for i in samples], file=fout)
            for file in expected_files:
                if file == "common":
                    cols = [true_str] * len(samples)
                else:
                    cols = [true_str if char == "1" else false_str for char in list(file)]
                print(file.ljust(wname, " "), *cols, file=fout)


rule bench:
    input:
        vcf = "3_filtered/{file}.final.vcf.gz",
        index = "3_filtered/{file}.final.vcf.gz.tbi"
    output:
        dir = directory("bench/{file}")
    log: "bench/{file}.log"
    params:
        sizefilt = int(minlength * 0.7),
        bed = f"--includebed {bench_bed}" if bench_bed is not None else ""
    shell:
        "truvari bench"
        " --base {bench_vcf}"
        " --comp {input.vcf}"
        " --output {output.dir}"
        " --sizefilt {params.sizefilt}"
        " --sizemin {minlength}"
        " --sizemax {maxlength}" 
        " {params.bed}"
        " {bench_args}"
        " 2> {log}"


rule vcf2db_bench:
    input:
        dir = "bench/{file}"
    output:
        joblib = "bench/{file}.jl"
    log: "bench/{file}.jl.log"
    shell:
        "truvari vcf2df -i -b {input.dir} {output.joblib} 2> {log}"


rule vcf2df_merge:
    input:
        vcf = "4_merged/merged.vcf"
    output:
        jl = "4_merged/merged.jl"
    log: "4_merged/merged.jl.log"
    params:
        samples_str = ",".join(samples)
    shell:
        "truvari vcf2df"
        " -i"
        " -m"
        " -f"
        " -s {params.samples_str}"
        " -S"
        " {input.vcf}"
        " {output.jl}"
        " 2> {log}"
