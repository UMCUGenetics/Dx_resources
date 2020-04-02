## Settings used in run_Exomdepth.py ##

# Location of repository
#cwd = "/hpc/diaggen/software/production/Dx_resources/ExomeDepth/"
cwd="/hpc/diaggen/software/development/Dx_resources_ED/ExomeDepth/"

# R script for RefSet and CNV calling
call_cnv_r = str(cwd)+"Exomedepth_callCNVs.R"
create_refset_r = str(cwd)+"Exomedepth_create_Ref_set.R"

#ExomeDepth CSV to VCF parser
csv2vcf = str(cwd)+"ed_csv_to_vcf.py"
vcf_template = str(cwd)+"template.vcf"

# Reference genome
reference_genome = "/hpc/diaggen/data/databases/ref_genomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"

#Used R version
r_version = "R/3.5.1"

# Location of target files
refset_dir = "/hpc/diaggen/data/databases/ExomeDepth_refset/"
#reffile_dir = "/hpc/diaggen/software/production/Dx_tracks/ExomeDepth/"
reffile_dir = "/hpc/diaggen/software/development/Dx_tracks_EDlfs/ExomeDepth/"

#Reference set
#refset = "Jan2020"
refset = "Nov2019"

analysis = {"HC":
              {"refset":{"female":str(refset_dir)+"/HC_female_"+str(refset)+".EDref",
                         "male":str(refset_dir)+"/HC_male_"+str(refset)+".EDref"},
               "target_bed":str(reffile_dir)+"High_confident_SureSelect_CREv2_elidS30409818_Covered_dp30_500_cv20.bed",
               "exon_bed":str(reffile_dir)+"exons.hg19.full_HC_CREv2_elidS30409818.tsv"
               },

          "UMCU": 
               {"refset":{"female":str(refset_dir)+"/UMCU_female_"+str(refset)+".EDref",
                          "male":str(refset_dir)+"/UMCU_male_"+str(refset)+".EDref"},
               "target_bed":str(reffile_dir)+"ENSEMBL_UCSC_merged_collapsed_sorted_v3_20bpflank_flat.bed",
               "exon_bed":str(reffile_dir)+"exons.hg19.full_UMCU20bp.tsv"
               }
         }   

# Parameters transition.probability
probability = {"HC": 0.0001, "UMCU": 0.5}

# qsub settings
qsub_call = "#!/bin/bash\n#$ -cwd\n#$ -pe threaded 1\n#$ -l h_vmem=10G\n#$ -l h_rt=6:00:00\n#$ -m a\n#$ -M "
qsub_ref = "#!/bin/bash\n#$ -cwd\n#$ -pe threaded 1\n#$ -l h_vmem=10G\n#$ -l h_rt=48:00:00\n#$ -m a\n#$ -M "

# Settings for GT, CN, and gender determination
gender = {"male":"M","female":"F"}
ratio_threshold_del = 0.25
par1 = [60001,2699520]
par2 = [154931044,155260560]
normal_CN = {"female":{"auto":2,"chrX":2,"chrXpar":2,"chrY":0},"male":{"auto":2,"chrX":1,"chrXpar":2,"chrY":1}}

# IGV session settings:
template_xml = str(cwd)+"igv_session_template.xml"
igv_xml = str(cwd)+"igv_xml_session.py"
igv_settings = {"ratio":[0, 1, 2], "log2ratio":[-2.5, 0, 2.5] }

