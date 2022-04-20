# Wrapper scripts for ExomeDepth (UMCU version)

All script are made tested using Python 3.6.8
## Installation 
To run the ExomeDepth wrapper scripts we need to create a virtual python environment and folder with all the required R packages.

### Make virtual python enviroment
```bash
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

### Build singularity R image 
Build a singularity container image from the [rocker tidyverse](https://hub.docker.com/layers/rocker/tidyverse/3.5.1/images/sha256-916d4e919fbac9ee1f06db7622b4731268dff0cb17998b3db6f629a3e58f73a7?context=explore) docker container.
```bash
singularity build rocker-tidyverse-3.5.1.img docker://rocker/tidyverse:3.5.1
```

### Install R packages in custom library location
Start a singularity shell using the [rocker tidyverse](https://hub.docker.com/layers/rocker/tidyverse/3.5.1/images/sha256-916d4e919fbac9ee1f06db7622b4731268dff0cb17998b3db6f629a3e58f73a7?context=explore) singularity image container, set R_LIBS location and open R. Make sure to change -B (root path for all required analysis files) and R_LIBS paths (location of R modules). 
```bash
singularity shell -B </change/this/path>:</change/this/path> /change/path/to/rocker-tidyverse-3.5.1.img
export R_LIBS=</change/path/to/R_libs/3.5.1>
R
```
Install R packages.
```R
library(devtools)
BiocManager::install(c("Biostrings","Rsamtools","GenomicRanges","VGAM","zlibbioc","bitops","BiocGenerics","BiocParallel","S4Vectors","IRanges","GenomeInfoDb","RCurl","GenomeInfoDbData","XVector","GenomicAlignments","Biobase","DelayedArray","matrixStats","Matrix","lattice","aod","stringr","stringi","dplyr","rlang","Rcpp","assertthat","glue","pkgconfig","tibble","pillar","crayon","vctrs","tidyselect","purrr","SummarizedExperiment","R6"))

install_github("UMCUGenetics/ExomeDepth")
```

## run_ExomeDepth.py
Main wrapper script for the ExomeDepth analysis.

First create a reference set which can be used to call CNV.

## Create refset
Refset will (default) create 4 output file in 4 seperate folders: a HighConfident (HC) and LowCondfident (UMCU) reference file (.EDref) for females and male seperate. Male and Female is determined on chrY read counts.
- Please do not include samples in the folder with known sex-chromsome deviations (e.g. XXY, X).
- The 4 .EDref files can be copied into a specific folder.
- Include this folder path and the reference set naming (i.e. Jan2020) in settings.py for CNV calling.
``` bash
. </path/to/repo>/venv/bin/activate
python run_ExomeDepth.py makeref <output folder> <inputfolder: folder containing realigned.BAM files> <prefix: i.e. Jan2020>
```

## Caling CNVs
CNV calling will result in a folder with a VCF file containing the significant CNV calls. For both the HC and UMCU target files seperately.

### Single sample
``` bash
. </path/to/repo>/venv/bin/activate
python run_ExomeDepth.py callcnv  <output folder> <input bam> <run id> <sample id> 
```

### Multisample (e.g. re-analysis of old run):
Use 2 threads for each sample. Thus to process 4 samples simultaneously use 8 threads.
``` bash
. </path/to/repo>/venv/bin/activate
python submit_batch_exomedepth.py <input folder> <output folder> <runid> <sample count>
```

## How to make a HC file
1) Start by making a stipped version of the original enrichment BED file:
	* remove headers
	* only include 3 columns: chromosome, start, stop
        * make sure no overlapping targets are present in the BED file (if present, merge overlapping fragments)
          This can be done with bedtools, for example: 
          ``` bash
          cat <original_bed> | sort -k 1V -k 2n -k 3n |   bedtools merge -i - > <flat_bed>
          ```
        * remove 'chr' from chromsome names, if needed

2) Optional: slice the stripped BED in fragments of (minimal) x bp. This is done to chop up large exons.
    ``` bash
    ./slice_bed_file.py <bed_file> <length> > <bed_file_sliced>
    ```

3) Make new HC BED file for two populations. Each consisting of minimal 50 males and 50 females.
    
    Make a folder for each population, and add BAM files into male and female specific folders, for example:
    Population1/male\
    Population1/female\
    Population2/male\
    Population2/female

    Calculate coverage stats for each male/female folder for each population
    ``` bash
    sh utils/run_sambamba.sh <full path to sambamba executable> <folder> <bed_file> <email>
    ```

    Calculate coverage stats overview for each male/female folder for each population
    ``` bash
    ./utils/filter_probe_file.py <folder> > <population>_<male/female>_output_all
    ```

4) Select least variable regions: 95% of chr1-22+chrX (female), and 33% of chrY (male)
    Calculate number of targets chr1-22+chrX
    ``` bash
    cat <bed_file> | awk '($1 != "Y")' | wc -l
    ```
    Calculate number of targets chrY
    ``` bash
    cat <bed_file> | awk '($1 == "Y")' | wc -l
    ```

5) Slice  95%/33% least variable regions
    ``` bash
    cat  <population>_female_output_all |sed 's/inf/99999/g'| awk '($1 != "Y")' | sort -nk6 | head -n <95% of chr1-22+chrX count> |sed 's/X/999999999/g' | sort -nk1 -nk2 |sed 's/999999999/X/g' | awk '{OFS="\t"; print $1,$2,$3,$4"_"$5"_"$6}' > <population>_female_auto_chrX
    cat  <population>_male_output_all |sed 's/inf/99999/g'| awk '($1 == "Y")' | sort -nk6 | head -n <33% of chrY count> | sort -nk1 -nk2 |awk '{OFS="\t"; print $1,$2,$3,$4"_"$5"_"$6}'> <population>_male_chrY
    cat <population>_female_auto_chrX <population>_male_chrY > <bed_file>_<population>
    ```

6) Calculate overlap bewteen the two populations. At the moment this should be >=99%.
    ``` bash
    cat <bed_file>_<population1> <bed_file>_<population2> | cut -f1,2,3 | sort | uniq -c | awk '($1==2)' |wc -l  # Overlap
    cat <bed_file>_<population1> | wc -l  # total targets (similar in both populations)
    ```
    Divide overlap/total targets. 
    Optional: if this is <99% step 4 can be repeated with decreasing the 95%, for example, 94%, etc. 

7) Make final BED and Exomedepth exon.tsv 
    ``` bash
    cat <bed_file>_<population1> <bed_file>_<population2> | cut -f1,2,3 | sort | uniq -c | awk '($1==2)' |  sed 's/ \+/\t/g'  |cut -f 3,4,5 | sed 's/X/999999999/g'| sed 's/Y/9999999999/g' | sort -nk1 -nk2 |sed 's/9999999999/Y/g' | sed 's/999999999/X/g' > <final_bed_file>
    cat <final_bed_file> |awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3 }' > exons.hg19.full.tsv
    ```
    Note: a tab-seperated header should be included in exons.hg19.full.tsv:
    ```
    chromosome	start	end	name
    ```
    Copy the final_bed_file and exons.hg19.full.tsv to a repository/location of choice.\
    Include in `setting.py` if these file are needed in the ExomeDepth analysis.

## Other scripts in this repository 
### check_gender_bam.py ###
Check gender of BAM files in a specific folder.\
Consider male is nonPAR regions of chrY fraction >=0.12 and chrX <= 2.3.\
Consider female is nonPAR regions of chrY fraction <=0.06 and chrX >= 4.0.\
All other combinations are considerend unknown (these include i.e. XXY and X or sample with mosaik/large CNV, and/or contamination).
``` bash
python check_gender_bam.py <inputfolder>
```

### identify_merge_samples.py
``` bash
python identify_merge_samples.py <inputfolder> <output name for list with merge_samples> 
```
This script will produce a list of merge samples based on slide barcode comparison in BAM and runID.\
Made specific for UMCU IAP and NF runs.

### make_BEDdetail.py
``` bash
python make_BEDdetail.py <inputfolder> <outputfolder> <output prefix> <list with merge_samples>
```
Make a BED detail file for UCSC browser. This BED file is a frequency summary of all events detected in the provided VCF files.\
WES-CNV files should be present in inputfolder.\
List with merge_samples is the output file of identify_merge_samples.py or tabseperated file inclusding sampleID and runID, in that order.

### exomedepth_summary.py 
Creates ExomeDepth summary stats file. This scripts need to be run with python3

### ed_csv_to_vcf.py
Converts ExomeDepth csv file to VCF using pyvcf

### filter_probe_file.py
Collected Sambamba coverage stats, and calculate mean, stdev, and CV.

### igv_xml_session.py
Makes IGV session xml files

### settings.py
This scripts contains path, files, and settings that are used in multiple python scripts

