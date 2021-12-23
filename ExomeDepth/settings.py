## Settings used in run_Exomedepth.py ##

# Location of repository
cwd = "/hpc/diaggen/software/production/Dx_resources/ExomeDepth"

# R script for RefSet and CNV calling
call_cnv_r = str(cwd)+"/Exomedepth_callCNVs.R"
create_refset_r = str(cwd)+"/Exomedepth_create_Ref_set.R"

#ExomeDepth CSV to VCF parser
csv2vcf = str(cwd)+"/ed_csv_to_vcf.py"
vcf_template = str(cwd)+"/template.vcf"

# HTML template
html = str(cwd)+"/template.html"

# Reference genome
reference_genome = "/hpc/diaggen/data/databases/ref_genomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"

#Used R version
r_library_path = "/hpc/diaggen/software/production/R_libs/ExomeDepthv2.0.1/3.5.1_singularity/"
singularity_r_container = "/hpc/diaggen/software/singularity_cache/rocker-tidyverse-3.5.1.img"
singularity_mount_path = "/hpc/diaggen:/hpc/diaggen"


# Location of target files
reffile_repo = "/hpc/diaggen/software/production/Dx_tracks"
reffile_dir = "{}/ExomeDepth".format(reffile_repo)

#Reference set
refset = "RS-SSv7-2021-1"		## Annotation for HC and UMCU refsets
callingmodel_HC = "HC_SSv7-2021-1"	## High Confident track annotation
callingmodel_UMCU = "UMCU"		## UMCU track annotation

analysis = {
    "HC":
        {"target_bed":"{0}/{1}_target.bed".format(reffile_dir, callingmodel_HC),
        "exon_bed":"{0}/{1}_exon.tsv".format(reffile_dir, callingmodel_HC),
        "calling_model":callingmodel_HC
        },

    "UMCU":
        {"target_bed":"{0}/{1}_target.bed".format(reffile_dir, callingmodel_UMCU),
        "exon_bed":"{0}/{1}_exon.tsv".format(reffile_dir, callingmodel_UMCU),
        "calling_model":callingmodel_UMCU
        }
    }


# Parameters transition.probability
probability = {"HC": 0.0001, "UMCU": 0.5}

# General settings
ratio_threshold_del = 0.25
par1 = [60001,2699520]
par2 = [154931044,155260560]
gender_determination_locus_y = 'Y:2649520-59034050'
gender_determination_locus_x = 'X:2699520-154931044'
gender_determination_y_ratio = [0.06, 0.12]
gender_determination_x_ratio = [2.3, 3.7]
normal_CN = {"female":{"auto":2,"chrX":2,"chrXpar":2,"chrY":0},"male":{"auto":2,"chrX":1,"chrXpar":2,"chrY":1}}
expectedCNVlength = 1000000
chromosome_order = {'1':0, '2':1, '3':2, '4':3, '5':4, '6':5, '7':6, '8':7, '9':8, '10':9, '11':10, '12':11, '13':12, '14':13, '15':14, '16':15, '17':16, '18':17, '19':18, '20':19, '21':20, '22':21, 'X':22, 'Y':23, 'MT':24}

# IGV session settings:
template_xml = str(cwd)+"/igv_session_template.xml"
igv_xml = str(cwd)+"/igv_xml_session.py"
igv_settings = {"ratio":[0, 1, 2], "log2ratio":[-2.5, 0, 2.5] }

# Check VCF stats criteria
correlation = 0.98
number_calls = [35,200]
del_dup_ratio = [15,85] 

reanalysis_dic = {
    "XXY":["female", "XXYMaleAsFemale"],
    "X": ["male", "XFemaleAsMale"],
    "asFemale": ["female", "AsFemale"],
    "asMale": ["male", "AsMale"]
    }

merge_warning = "DO_NOT_USE_MergeSample"

