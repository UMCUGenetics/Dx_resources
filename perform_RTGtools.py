#! /usr/bin/env python
import sys, os, re
import commands
import getopt
from optparse import OptionParser
from optparse import OptionGroup
import random,string

## Author: M.G. Elferink
## Date: 05-10-2016
## Purpose: Perform RTGtools analyses for VKLG-GIAB samples.
## Usage: ./perform_RTGtools.py -i [vcf] 

def SNP_s():
	sh.write("bcftools filter "+str(vcf)+" -i \'TYPE=\"snp\"\'| sed 's/^chr//g' > "+str(vcf)[0:-4]+"_SNP.vcf\n")  # And remove "chr"	
	sh.write("bgzip "+str(vcf)[0:-4]+"_SNP.vcf\n")
	sh.write("tabix -f "+str(vcf)[0:-4]+"_SNP.vcf.gz\n")

def INDEL_s():
	sh.write("bcftools filter "+str(vcf)+" -i \'TYPE=\"indel\"\'| sed 's/^chr//g' > "+str(vcf)[0:-4]+"_INDEL.vcf\n") # And remove "chr"
	sh.write("bgzip "+str(vcf)[0:-4]+"_INDEL.vcf\n")
     	sh.write("tabix -f "+str(vcf)[0:-4]+"_INDEL.vcf.gz\n")

def SNP_i():
	sh.write("bedtools intersect -a "+str(vcf)[0:-4]+"_SNP.vcf.gz -b "+str(high_conf_bed) +" -u -header > " +str(vcf)[0:-4]+"_SNP_intersect.vcf\n")
	sh.write("bgzip "+str(vcf)[0:-4]+"_SNP_intersect.vcf\n")
	sh.write("tabix -f "+str(vcf)[0:-4]+"_SNP_intersect.vcf.gz\n")

def INDEL_i():
	sh.write("bedtools intersect -a "+str(vcf)[0:-4]+"_INDEL.vcf.gz -b "+str(high_conf_bed) +" -u -header > " +str(vcf)[0:-4]+"_INDEL_intersect.vcf\n")
	sh.write("bgzip "+str(vcf)[0:-4]+"_INDEL_intersect.vcf\n")
	sh.write("tabix -f "+str(vcf)[0:-4]+"_INDEL_intersect.vcf.gz\n")

def SNP_r():
	## All records
	sh.write("echo RTG calculations SNP ALL variants\n")
	sh.write("java -jar -Xmx10G "+str(rtg)+" vcfeval -t " + str(SDF) + " -T 6 --sample="  +str(sample_giab)+","+str(sample_test)+" --baseline="+str(giab_snp)+" --calls="+str(vcf)[0:-4]+"_SNP_intersect.vcf.gz "+ " --output="+str(vcf)+"_SNP_intersect_ALL"+" --bed-regions="+str(bed)+" --all-record\n")
	## PASS only
	sh.write("echo RTG calculations SNP PASS variants\n")
	sh.write("java -jar -Xmx10G "+str(rtg)+" vcfeval -t " + str(SDF) + " -T 6 --sample="  +str(sample_giab)+","+str(sample_test)+" --baseline="+str(giab_snp)+" --calls="+str(vcf)[0:-4]+"_SNP_intersect.vcf.gz "+ " --output="+str(vcf)+"_SNP_intersect_PASS"+" --bed-regions="+str(bed)+"\n")

def INDEL_r():
	## All records
	sh.write("echo RTG calculations INDEL ALL variants\n")
	sh.write("java -jar -Xmx10G "+str(rtg)+" vcfeval -t " + str(SDF) + " -T 6 --sample="  +str(sample_giab)+","+str(sample_test)+" --baseline="+str(giab_indel)+" --calls="+str(vcf)[0:-4]+"_INDEL_intersect.vcf.gz "+ " --output="+str(vcf)+"_INDEL_intersect_ALL"+" --bed-regions="+str(bed)+" --all-record\n")
	## PASS only	
	sh.write("echo RTG calculations INDEL PASS variants\n")	
	sh.write("java -jar -Xmx10G "+str(rtg)+" vcfeval -t " + str(SDF) + " -T 6 --sample="  +str(sample_giab)+","+str(sample_test)+" --baseline="+str(giab_indel)+" --calls="+str(vcf)[0:-4]+"_INDEL_intersect.vcf.gz "+ " --output="+str(vcf)+"_INDEL_intersect_PASS"+" --bed-regions="+str(bed)+"\n")

################
if __name__ == "__main__":
        parser = OptionParser();
        group = OptionGroup(parser, "Main options")
	group.add_option("-i", dest="vcf", type='string', help="VCF input file")
	group.add_option("-o", default="combined", dest="option",choices=['combined', 'indel', 'snp'], help="options for analysis [combined VCF [\"combined\", default], indel VCF [indel], snp-only VCF [snp]]")
	group.add_option("-r", default="/hpc/local/CentOS7/cog_bioinf/rtg-tools-3.6.2/RTG.jar", dest="RTG", help="RTG path [default /hpc/local/CentOS7/cog_bioinf/rtg-tools-3.6.2//RTG.jar]")
	group.add_option("-c", default="/hpc/cog_bioinf/common_dbs/GIAB/NIST_2.19/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio_noMT.bed", dest="high_conf_bed", help="High confident BED GIAB file [default /hpc/cog_bioinf/common_dbs/GIAB/NIST_2.19/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio_noMT.bed]")
	group.add_option("-b", default="/hpc/cog_bioinf/diagnostiek/production/Dx_resources/Tracks/ENSEMBL_UCSC_merged_collapsed_sorted_v2_20bpflank.bed", dest="target_bed", help="Target BED file [default /hpc/cog_bioinf/diagnostiek/production/Dx_resources/Tracks/ENSEMBL_UCSC_merged_collapsed_sorted_v2_20bpflank.bed ]")
	group.add_option("--sdf", default="/hpc/local/CentOS7/cog_bioinf/rtg-tools-3.6.2/Homo_sapiens.GRCh37.GATK.illumina.SDF", dest="SDF", help="SDF file [default /hpc/local/CentOS7/cog_bioinf/rtg-tools-3.6.2/Homo_sapiens.GRCh37.GATK.illumina.SDF ]")
	group.add_option("--snp", default="/hpc/cog_bioinf/common_dbs/GIAB/NIST_2.19/SNP_GIAB12878_nist2.19_truth.vcf.gz", dest="giab_snp", help="GIAB SNP truth VCF [default /hpc/cog_bioinf/common_dbs/GIAB/NIST_2.19/SNP_GIAB12878_nist2.19_truth.vcf.gz]")
	group.add_option("--indel", default="/hpc/cog_bioinf/common_dbs/GIAB/NIST_2.19/INDELS_GIAB12878_nist2.19_truth.vcf.gz", dest="giab_indel", help="GIAB INDEL truth VCF [default /hpc/cog_bioinf/common_dbs/GIAB/NIST_2.19/INDEL_GIAB12878_nist2.19_truth.vcf.gz]")
	group.add_option("-t", default="GIAB12878", dest="sample_test", help="sample name for Test sample [default GIAB12878]")
	group.add_option("-g", default="GIAB12878", dest="sample_giab", help="Target BED file GIAB sample [default GIAB12878]")
	parser.add_option_group(group)
        (opt, args) = parser.parse_args()

	if opt.vcf:
		vcf=str(opt.vcf)
	else:	
		sys.exit("provide VCF [-i VCF ]")

	if (os.path.exists(str(vcf))):
		pass
	else:
		sys.exit("VCF does not exist")

	option=str(opt.option)
	if option == "combined" or option == "indel" or option == "snp":
		pass
	else:
		sys.exit("Option not correct [-o [combined/indel/snp]")

	print "VCF file = "+str(vcf)+"\nOption = "+str(option)
	
	high_conf_bed=str(opt.high_conf_bed)
	rtg=str(opt.RTG)
	bed=str(opt.target_bed)
	SDF=str(opt.SDF)
	sample_test=str(opt.sample_test)
	sample_giab=str(opt.sample_giab)
	giab_snp=str(opt.giab_snp)
	giab_indel=str(opt.giab_indel)

	tag=''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6))
	sh=open("Submit_RTG_"+str(tag)+".sh","w")
	sh.write("module load bed/bedtools/2.25.0\n")
	sh.write("module load tabix/0.2.6\n")
	sh.write("module load vcfbcf/bcftools/1.3\n")

	## Split VCF in INDELs and SNPs
	if option == "combined":
		SNP_s()
		INDEL_s()
	if option == "indel": 	## just to make sure all VCF files are treated equal!
		INDEL_s()
	if option =="snp": 	## just to make sure all VCF files are treated equal!
		SNP_s()

	## Intersect VCF files with high confident BED file
	if option == "combined":
		SNP_i()
		INDEL_i()
	if option == "indel": 
		INDEL_i()
	if option =="snp": 
		SNP_i()

	## Perform RTGtools comparison
	if option == "combined":
	        SNP_r()
	        INDEL_r()
	if option == "indel":
	        INDEL_r()
	if option =="snp":
		SNP_r()
	sh.close()

	## sumbit to cluster ##
	print "qsub -cwd -q all.q -l h_rt=01:00:00 -l h_vmem=20G -pe threaded 6 Submit_RTG_"+str(tag)+".sh"
	os.system("qsub -cwd -q all.q -l h_rt=01:00:00 -l h_vmem=20G -pe threaded 6 Submit_RTG_"+str(tag)+".sh" )
	
	print "finished"
