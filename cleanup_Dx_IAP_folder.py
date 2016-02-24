#! /usr/bin/env python
import sys, os, re
import commands
import getopt

## Log git commits of used repo's
os.system(" date >> GIT_IAP.log")
os.system(" date >> GIT_ExonCov.log")
os.system(" date >> GIT_Dx_resources.log")
os.system("git --git-dir=/hpc/cog_bioinf/common_scripts/IAP_Dx_TARGET/.git log > GIT_IAP.log")
os.system("git --git-dir=/hpc/cog_bioinf/common_scripts/ExonCov/.git log > GIT_ExonCov.log")
os.system("git --git-dir=/hpc/cog_bioinf/data/mapping/diagnostiek/Dx_resources/.git log > GIT_Dx_resources.log")

# Make run_stats.txt file from all flagstats
os.system("/hpc/cog_bioinf/data/mapping/diagnostiek/Dx_resources/get_stats_from_flagstat.pl >run_stats.txt")

# move unused VCFs in Redundant folder
os.system("mkdir Redundant_VCF_files/")
os.system("find -iname \"*phased.vcf*\" -exec mv {} Redundant_VCF_files/ \;")
os.system("find -iname \"*raw_variants.vcf*\" -exec mv {} Redundant_VCF_files/ \;")

# cp .bai into bam.bai
pwd= commands.getoutput("pwd")
action="find "+str(pwd)+" -iname \"*dedup.realigned.bai\"" 
files=commands.getoutput(action).split()
for file in files:
	action= "cp " + file+ " "+ file[:-1]+"m.bai"
	print action
	os.system(action)


sys.exit("Finished")

