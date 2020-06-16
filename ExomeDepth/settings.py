## Settings used in run_Exomdepth.py ##

# Location of repository
cwd = "/hpc/diaggen/software/production/Dx_resources/ExomeDepth/"

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
reffile_dir = "/hpc/diaggen/software/production/Dx_tracks/ExomeDepth/"

#Reference set
refset = "May2020v2"

analysis = {
    "HC":
        {"target_bed":str(reffile_dir)+"High_confident_SureSelect_CREv2_elidS30409818_Covered_300bp_95_overlap_v2.bed",
        "exon_bed":str(reffile_dir)+"exons.hg19.full_300bp_95_overlap_v2.tsv"
        },

    "UMCU":
        {"target_bed":str(reffile_dir)+"ENSEMBL_UCSC_merged_collapsed_sorted_v3_20bpflank_flat_v1.bed",
        "exon_bed":str(reffile_dir)+"exons.hg19.full_UMCU20bp_v1.tsv"
        }
    }


# Parameters transition.probability
probability = {"HC": 0.0001, "UMCU": 0.5}

# General settings
ratio_threshold_del = 0.25
par1 = [60001,2699520]
par2 = [154931044,155260560]
normal_CN = {"female":{"auto":2,"chrX":2,"chrXpar":2,"chrY":0},"male":{"auto":2,"chrX":1,"chrXpar":2,"chrY":1}}
expectedCNVlength = 1000000

# IGV session settings:
template_xml = str(cwd)+"igv_session_template.xml"
igv_xml = str(cwd)+"igv_xml_session.py"
igv_settings = {"ratio":[0, 1, 2], "log2ratio":[-2.5, 0, 2.5] }

