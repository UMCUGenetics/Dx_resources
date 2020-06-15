# Wrapper scripts for ExomeDepth (UMCU version)

All script are made tested using Python 3.6.8
## Make virtual python enviroment
``` bash
virtualenv -p python3 venv
. venv/bin/activate
easy_install pip
pip install -r requirements.txt
```

## run_ExomeDepth.py
Main wrapper script for the ExomeDepth analysis\

First create a reference set which can be used to call CNV

## Create refset
``` bash
. {path_to_repo}/venv/bin/activate
python run_ExomeDepth.py makeref {output folder} {inputfolder: folder containing realigned.BAM files} {prefix: i.e. Jan2020}
```
Refset will (default) create 4 output file in 4 seperate folders: a HighConfident (HC) and LowCondfident (UMCU) reference file (.EDref) for females and male seperate. \
Male and Female is determined on chrY read counts.\
Please do not include samples in the folder with known sex-chromsome deviations (e.g. XXY, X)\
The 4 .EDref files can be copied into a specific folder.\
Include this folder path + the reference set naming (i.e. Jan2020) in settings.py for CNV calling.

## Caling CNVs
Single sample
``` bash
. {path_to_repo}/venv/bin/activate
python run_ExomeDepth.py callcnv  {output folder} {input bam} {run id} {sample id} {refset prefix: i.e. Jan2020)
```
CNV calling will result in a folder with a VCF file containing the significant CNV calls. For both the HC and UMCU target files seperately

Multisample (e.g. re-analysis of old run):\
NOTE: Use 2 threads for each sample. Thus to process 4 samples simulatnious use 8 threads

``` bash
. {path_to_repo}/venv/bin/activate
python submit_batch_exomedepth.py {input folder} {output folder} {samples(/threads)}
```




## How to make a HC file
1) Optional: slice the original BED in fragments of x bp. This is done to chop up larger exons
``` bash
./slice_bed_file.py {bed_file} {length} > {bed_file_sliced}
```

2) Make new HC BED file for two populations. Each consisting of minimal 50 males and 50 females.
Make a folder for each population, and add BAM files into male and female folders.

Calculate coverage stats for each male/female folder for each population
``` bash
sh run_sambamba.sh {folder} {bed_file} {email}
```
Calculate coverage stats overview for each male/female folder for each population
``` bash
./filter_probe_file_all.py {folder} > {population}_{male/female}_output_all
```

3) Select least variable regions: 95% of chr1-22+chrX (female), and 33% of chrY (male)
Calculate number of targets chr1-22+chrX
``` bash
cat {bed_file} | awk '($1 != "Y")' | wc -l
```
Calculate number of targets chrY
``` bash
cat {bed_file} | awk '($1 == "Y")' | wc -l
```

4) Slice  95%/33% least variable regions
``` bash
cat  {population}_female_output_all |sed 's/inf/99999/g'| awk '($1 != "Y")' | sort -nk6 | head -n {95% of chr1-22+chrX count} |sed 's/X/999999999/g' | sort -nk1 -nk2 |sed 's/999999999/X/g' | awk '{OFS="\t"; print $1,$2,$3,$4"_"$5"_"$6}' > {population}_female_auto_chrX
cat  {population}_male_output_all |sed 's/inf/99999/g'| awk '($1 == "Y")' | sort -nk6 | head -n {33% of chrY count} | sort -nk1 -nk2 |awk '{OFS="\t"; print $1,$2,$3,$4"_"$5"_"$6}'> {population}_male_chrY
cat {population}_female_auto_chrX {population}_male_chrY > {bed_file}_{population}
```

5) Calculate overlap bewteen the two populations. This should be minimal 99%.
``` bash
cat {bed_file}_{population1} {bed_file}_{population2} | cut -f1,2,3 | sort | uniq -c | awk '($1==2)' |wc -l  # Overlap
cat {bed_file}_{population1} | wc -l  # total targets (similar in both populations)
```
Divide overlap/total targets. 
Optional: if this is <99% step 3 can be repeated with decreasing the 95%, for example, 94%, etc.


6) Make final BED and Exomedepth exon.tsv 
``` bash
cat {bed_file}_{population1} {bed_file}_{population2} | cut -f1,2,3 | sort | uniq -c | awk '($1==2)' |  sed 's/ \+/\t/g'  |cut -f 3,4,5 | sed 's/X/999999999/g'| sed 's/Y/9999999999/g' | sort -nk1 -nk2 |sed 's/9999999999/Y/g' | sed 's/999999999/X/g' > {final_bed_file}
cat {final_bed_file} |awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3 }' > exons.hg19.full.tsv
```

Copy the final_bed_file and exons.hg19.full.tsv to a repository/location of choice.\
Include in setting.py if these file are needed in the ExomeDepth analysis.\



## Other scripts in this repository 
#### ed_csv_to_vcf.py
Converts ExomeDepth csv file to VCF using pyvcf

#### igv_xml_session.py
Makes IGV session xml files

#### settings.py
This scripts contains path, files, and settings that are used in multiple python scripts

