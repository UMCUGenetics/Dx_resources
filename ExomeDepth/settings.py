## Settings used in run_Exomdepth.py ##

# Location of repository
cwd = "/hpc/cog_bioinf/diagnostiek/development/Dx_resources_ED/ExomeDepth/"
# R script for RefSet and CNV calling
call_CNV=str(cwd)+"Exomedepth_callCNVs.R"
refSet=str(cwd)+"Exomedepth_create_Ref_set.R"
#ExomeDepth CSV to VCF parser
csv2vcf= str(cwd)+"ed_csv_to_vcf.py"
template=str(cwd)+"template.vcf"

# Reference genome
reference_genome= "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"

#Used R version
Rversion="R/3.5.1"

# Location of target files
refset_dir="/hpc/cog_bioinf/diagnostiek/projects/WES_CNV_ExomeDepth/Double_reference_analyses/RefSet/Final_RefSet/Full_set/"

#Reference set
refset="Oct2019"
analysis={"HC":
              {"refset":{"female":str(refset_dir)+"/HC_female_"+str(refset)+".EDref",
                         "male":str(refset_dir)+"/HC_male_"+str(refset)+".EDref"},
               "target_bed":str(refset_dir)+"High_confident_SureSelect_CREv2_elidS30409818_Covered_dp30_500_cv20.bed",
               "exon_bed":str(refset_dir)+"exons.hg19.full_HC_CREv2_elidS30409818.tsv"
               },

          "UMCU": 
               {"refset":{"female":str(refset_dir)+"/UMCU_female_"+str(refset)+".EDref",
                          "male":str(refset_dir)+"/UMCU_male_"+str(refset)+".EDref"},
               "target_bed":str(refset_dir)+"ENSEMBL_UCSC_merged_collapsed_sorted_v3_20bpflank_flat.bed",
               "exon_bed":str(refset_dir)+"exons.hg19.full_UMCU20bp.tsv"
               }
         }   

# Parameters transition.probability
probability={"HC": 0.0001, "UMCU": 0.5}

# qsub settings
qsub_call="#!/bin/bash\n#$ -cwd\n#$ -pe threaded 1\n#$ -l h_vmem=10G\n#$ -l h_rt=6:00:00\n#$ -m a\n#$ -M "
qsub_ref="#!/bin/bash\n#$ -cwd\n#$ -pe threaded 1\n#$ -l h_vmem=10G\n#$ -l h_rt=48:00:00\n#$ -m a\n#$ -M "

# Settings for GT, CN, and gender determination
gender={"male":"M","female":"F"}
ratio_threshold_del=0.25
PAR1=[10001,2781479]
PAR2=[155701383,156030895]
normal_CN={"female":{"auto":2,"chrX":2,"chrXPAR":2,"chrY":0},"male":{"auto":2,"chrX":1,"chrXPAR":2,"chrY":1}}

