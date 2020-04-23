# Wrapper scripts for ExomeDepth (UMCU version)

All script are made tested using Python 3.6.8
## Make virtual python enviroment
virtualenv -p python3 venv
. venv/bin/activate
easy_install pip
pip install -r requirements.txt


## run_ExomeDepth.py
Main wrapper script for the ExomeDepth analysis\
run_ExomeDepth.py uses SGE (UMCU HPC enviroment) for reference set creation and parallel multisample CNV-calling.\
In addition there is an option for serial single sample CNV-calling that does not use SGE. Setial single sample processing is not avaiable for reference set creation.

First create a reference set which can be used to call CNV

## Create refset
``` bash
python run_ExomeDepth.py makeref {output folder} {inputfolder: folder containing realigned.BAM files} {prefix: i.e. Jan2020}
```
Refset will (default) create 4 output file in 4 seperate folders: a HighConfident (HC) and LowCondfident (UMCU) reference file (.EDref) for females and male seperate. \
Male and Female is determined on chrY read counts.\
Please do not include samples in the folder with known sex-chromsome deviations (e.g. XXY, X)\
The 4 .EDref files can be copied into a specific folder.\
Include this folder, and the reference set naming in settings.py for CNV calling.

## Caling CNVs
Single sample
``` bash
python run_ExomeDepth.py callcnv  {output folder} {input bam} {run id} {sample id} {refset prefix: i.e. Jan2020)
```
CNV calling will result in a folder with a VCF file containing the significant CNV calls. For both the HC and UMCU target files seperately

Multisample (e.g. re-analysis of old run):
NOTE: use the same amount of threads as samples!
NOTE: python3 enviroment is needed

``` bash
. {path_to_repo}/venv/bin/activate
python submit_batch_exomedepth.py {input folder} {output folder} {samples(/threads)}
```

## How to make a HC file
Make seperate folders for male and female, and sotflink/copy minimum of 50 BAM files (each) into these folders.\
  Calculate coverage statistics for target BED file (e.g. SureSelect_CREv2_elidS30409818_Covered.bed) of each BAM using Sambamba:
  ``` bash
  sh run_sambamba.sh {input_folder} {enrichment bed file}
  ```
  This step will result in a .bam_coverage file for each BAM.

  Next, filter variable, low or high coverage region in both the female and male folders:
  ``` bash
  python filter_probe_file.py {input_folder} {min coverage} {max coverage} {max coefficient of variation} > {output_filtering}
  ```
  e.g.
  ``` bash
       python ./filter_probe_file.py ./males/ 30 500 20 > male_output_30-500_20
       python ./filter_probe_file.py ./females/ 30 500 20 > female_output_30-500_20
  ```

  Calculate overlapping regions in male and female:

  Make autosomal High_confident file:
  ``` bash
  cat {output_filtering_female} {output_filtering_male} | cut -f1,2,3 |sort |uniq -c | awk '($1==2)'| sed 's/ /\t/g' | sed 's/\t\t/\t/g' |cut -f5,6,7 | awk '($1 != "X" && $1 != "Y")' |sort -nk1 -nk2 > {output_autosomal}
  ```
  e.g.
  ``` bash
  cat female_output_30-500_20 male_output_30-500_20| cut -f1,2,3 |sort |uniq -c | awk '($1==2)'| sed 's/ /\t/g' | sed 's/\t\t/\t/g' |cut -f5,6,7 | awk '($1 != "X" && $1 != "Y")' |sort -nk1 -nk2 > High_confident_SureSelect_CREv2_elidS30409818_Covered_dp30_500_cv20_noSex.bed
  ```

  Make chrX and chrY high confident based on males only:
  ``` bash
  cat {output_filtering_male} | sed 's/ //g' | awk '( $1== "X" || $1=="Y")' |cut -f1,2,3 > {output_sexchromosome}
  ```
  e.g.
  ``` bash
  cat male_output_30-500_20 | sed 's/ //g' | awk '( $1== "X" || $1=="Y")' |cut -f1,2,3 > High_confident_SureSelect_CREv2_elidS30409818_Covered_dp30_500_cv20_male.bed
  ```
  Note: threshold for coverage are similar to autosomal which might result in overfiltering on chrX and chrY becaus males are cn = 1.
 
  Make final High confident bed file:
  ``` bash
  cat {output_autosomal} {output_sexchromosome} > {output_final_bed}
  ```
  e.g.
  ``` bash
  cat High_confident_SureSelect_CREv2_elidS30409818_Covered_dp30_500_cv20_noSex.bed High_confident_SureSelect_CREv2_elidS30409818_Covered_dp30_500_cv20_male.bed > High_confident_SureSelect_CREv2_elidS30409818_Covered_dp30_500_cv20.bed
  ``` 

  Finally, as Exomdepth requires a slightly different format for targets, we need to make an 'exon.hg19' file:
  ``` bash
  echo -e "chromosome\tstart\tend\tname" > header
  cat {output_final_bed} | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3}' > {output_final_bed_4kol}
  cat header {output_final_bed_4kol} > {exon.hg19_file}
  ```
  e.g.
  ``` bash
  echo -e "chromosome\tstart\tend\tname" > header
  cat High_confident_SureSelect_CREv2_elidS30409818_Covered_dp30_500_cv20.bed | awk '{OFS="\t"; print $1,$2,$3,$1":"$2"-"$3}' > High_confident_SureSelect_CREv2_elidS30409818_Covered_dp30_500_cv20_4kol.bed
  cat header High_confident_SureSelect_CREv2_elidS30409818_Covered_dp30_500_cv20_4kol.bed > exons.hg19.full_HC_CREv2_elidS30409818.tsv
  ```

  Copy the HC bed file {output_final_bed} and exon.h19 {exon.hg19_file} to a repository/location of choice.\
  Include in setting.py if these file are needed in the ExomeDepth analysis.\

## Other scripts in this repository 
#### ed_csv_to_vcf.py
Converts ExomeDepth csv file to VCF using pyvcf

#### igv_xml_session.py
Makes IGV session xml files

#### settings.py
This scripts contains path, files, and settings that are used in multiple python scripts

