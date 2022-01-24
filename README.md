# lsvtool

A small tool to perform comparisons between large structural variants (SVs) from different callers and/or datasets.
The input is VCF (`*.vcf`/`*.vcf.gz`) files which will be normalized, filtered before being intersected to generate
figures.

**Tools used:**

- [`SURVIVOR`](https://github.com/fritzsedlazeck/SURVIVOR)
- [`truvari`](https://github.com/spiralgenetics/truvari)
- [`bcftools`](https://samtools.github.io/bcftools/bcftools.html)
- [`bedtools`](https://bedtools.readthedocs.io/en/latest/index.html)

**Steps:**

- Select non-filtered SV calls (labeled `PASS`) of correct type (e.g. `DEL`, `INS`, `INV` ...)
- Collapse proximal SV calls using `truvari collapse` to avoid superfluous calls
- Filter SVs against supplied blacklist, segmental duplications and for size.
- Compare and merge SVs between samples to find intersections.
- Generate heatmaps and venn diagram over sample intersects.
- (OPTIONAL) Benchmark calls against benchmark set of SVs.
- (OPTIONAL) Generate batch files for `igv` covering regions with the different unique and shared SVs for the samples. Used to plot support in covered regions.

## Installation and Running

The installation requires conda environment which will collect all the needed software. So, first make sure you have a working conda!

### Installation

    git clone https://github.com/HSiga/lsvtool.git  
    cd lsvtool  
    conda env create -n lsv --file env.yaml 
    conda activate lsv  
    pip install .  
    lsvtool -h

### Running analysis

Initialize a working directory

    lsvtool init -i file1.vcf.gz file2.vcf.gz.. -o output_directory  
    cd output_directory 

Adjust the parameters in the config (`parameters.config`) file. Finally to run the analysis use_

    lsvtool run
