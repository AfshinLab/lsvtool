# lsvtool
A small tool to perform comparisons between different large structural variants files generated by several structural variants callers.

The input files are <code>vcf.gz</code> files which will be interesected with each other.

The output will be: <br/>
    - Filtered <code>vcf/bedpe</code> files based on the different cutoffs. <br/>
    - A merged <code>vcf/bedpe</code> file with a comparison matrix, and a figure of intersecting regions between the inputs. <br/>
    - Lists of bedpe files for common and uniq regions for the filtered inputs, which can be easily plotted on IGV.
    Addtionally, Batch files which can be used to blot the regions in IGV after loading the bam files of the reads.


# Installation and Running
The installation requires conda environment which will collect all the needed software. So, first make sure you have a working conda!
## Installation:
<code>git clone https://github.com/HSiga/lsvtool.git </code><br/> 
<code>cd lsvtool  </code><br/>
<code>conda env create -n lsv2 --file env.yaml </code><br/>
<code>conda activate lsv2  </code><br/>
<code>pip install .  </code><br/>
<code>lsvtool -h</code>

## Running analysis:
<code>lsvtool init -i file1.vcf.gz file2.vcf.gz.. -o output_directory </code><br/> 
<code>cd output_directory </code>

Then adjust the parameters in the config file, and finally.

<code>lsvtool run</code>

