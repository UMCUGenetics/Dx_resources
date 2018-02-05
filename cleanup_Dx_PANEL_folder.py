#! /usr/bin/env python
import sys, os, re
import commands
import getopt

# Input folders
IAP_folder="/hpc/local/CentOS7/cog_bioinf/IAP_Dx/"
ExonCov_folder="/hpc/local/CentOS7/cog_bioinf/ExonCov/"
Dx_resources_folder="/hpc/cog_bioinf/diagnostiek/production/Dx_resources/"
Dx_INI_folder="/hpc/cog_bioinf/diagnostiek/production/Dx_INI/"
Dx_tracks_folder="/hpc/cog_bioinf/diagnostiek/production/Dx_tracks/"
# Get GIT tag version
IAP_v=commands.getoutput("git --git-dir="+str(IAP_folder)+".git describe --tags")
ExonCov_v=commands.getoutput("git --git-dir="+str(ExonCov_folder)+".git describe --tags")
Dx_resources_v=commands.getoutput("git --git-dir="+str(Dx_resources_folder)+".git describe --tags")
Dx_INI_v=commands.getoutput("git --git-dir="+str(Dx_INI_folder)+".git describe --tags")
Dx_tracks_v=commands.getoutput("git --git-dir="+str(Dx_tracks_folder)+".git describe --tags")
# Output files
output_IAP="GIT_IAP_"+str(IAP_v)+".log"
output_ExonCov="GIT_ExonCov_"+str(ExonCov_v)+".log"
output_Dx_resources="GIT_Dx_resources_"+str(Dx_resources_v)+".log"
output_Dx_INI="GIT_Dx_INI_"+str(Dx_INI_v)+".log"
output_Dx_tracks="GIT_Dx_tracks_"+str(Dx_tracks_v)+".log"

## Log git commits of used repo's
os.system(" date >> "+str(output_IAP))
os.system(" date >> "+str(output_ExonCov))
os.system(" date >> "+str(output_Dx_resources))
os.system(" date >> "+str(output_Dx_INI))
os.system(" date >> "+str(output_Dx_tracks))
## Print used folder in log files
os.system("echo "+str(IAP_folder)+" >> "+str(output_IAP))
os.system("echo "+str(ExonCov_folder)+" >> "+str(output_ExonCov))
os.system("echo "+str(Dx_resources_folder)+" >> "+str(output_Dx_resources))
os.system("echo "+str(Dx_INI_folder)+" >> "+str(output_Dx_INI))
os.system("echo "+str(Dx_tracks_folder)+" >> "+str(output_Dx_tracks))
## Print GIT version in log files
os.system("echo "+str(IAP_v)+" >> "+str(output_IAP))
os.system("echo "+str(ExonCov_v)+" >> "+str(output_ExonCov))
os.system("echo "+str(Dx_resources_v)+" >> "+str(output_Dx_resources))
os.system("echo "+str(Dx_INI_v)+" >> "+str(output_Dx_INI))
os.system("echo "+str(Dx_tracks_v)+" >> "+str(output_Dx_tracks))
## Print full GIT log in log files
os.system("git --git-dir="+str(IAP_folder)+".git log >> "+str(output_IAP))
os.system("git --git-dir="+str(ExonCov_folder)+".git log >> "+str(output_ExonCov))
os.system("git --git-dir="+str(Dx_resources_folder)+".git log >> "+str(output_Dx_resources))
os.system("git --git-dir="+str(Dx_INI_folder)+".git log >> "+str(output_Dx_INI))
os.system("git --git-dir="+str(Dx_tracks_folder)+".git log >> "+str(output_Dx_tracks))

# Make run_stats.txt file from all flagstats
os.system("/hpc/cog_bioinf/diagnostiek/production/Dx_resources/get_stats_from_flagstat.pl >run_stats.txt")

# Upload run data to trend analysis database
print "Uploading run data to trend analysis database\n"
pwd = commands.getoutput("pwd")
trend_analysis_command = ". {trend_analysis_path}/venv/bin/activate && python {trend_analysis_path}/trend_analysis.py upload processed_data {run_folder}".format(
    trend_analysis_path=' /hpc/cog_bioinf/diagnostiek/development/Trend_Analysis_tool',
    run_folder=pwd
)
os.system(trend_analysis_command)

# move unused VCFs in Redundant folder
os.system("mkdir Redundant_VCF_files/")
#os.system("find -iname \"*phased.vcf*\" -exec mv {} Redundant_VCF_files/ \;")
os.system("find -iname \"*raw_variants.vcf*\" -exec mv {} Redundant_VCF_files/ \;")

# cp .bai into bam.bai
pwd= commands.getoutput("pwd")
action="find "+str(pwd)+" -iname \"*dedup.realigned.bai\"" 
files=commands.getoutput(action).split()
for file in files:
	action= "cp " + file+ " "+ file[:-1]+"m.bai"
	print action
	os.system(action)

# change dots in Exoncov files into comma's
action2=" sed -i 's/\./\,/g' "+str(pwd)+"/Exoncov*/*tsv" 	## include use Exoncov output folder?
action3=" sed -i 's/\./\,/g' "+str(pwd)+"/Exoncov*/*/*tsv"	## include use Exoncov output folder?
os.system(action2)
os.system(action3)

sys.exit("Finished")

