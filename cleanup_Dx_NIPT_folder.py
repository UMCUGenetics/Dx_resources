#! /usr/bin/env python
import sys, os, re
import commands
import getopt

IAP_folder="/hpc/local/CentOS7/cog_bioinf/IAP_Dx/"
Chromate_folder="/hpc/local/CentOS7/cog_bioinf/chromate/"
Dx_resources_folder="/hpc/cog_bioinf/diagnostiek/production/Dx_resources/"
Dx_INI_folder="/hpc/cog_bioinf/diagnostiek/production/Dx_INI/"
# Get GIT tag version
IAP_v=commands.getoutput("git --git-dir="+str(IAP_folder)+".git describe --tags")
Chromate_v=commands.getoutput("git --git-dir="+str(Chromate_folder)+".git describe --tags")
Dx_resources_v=commands.getoutput("git --git-dir="+str(Dx_resources_folder)+".git describe --tags")
Dx_INI_v=commands.getoutput("git --git-dir="+str(Dx_INI_folder)+".git describe --tags")
# Output files
output_IAP="GIT_IAP_"+str(IAP_v)+".log"
output_Chromate="GIT_chromate_"+str(Chromate_v)+".log"
output_Dx_resources="GIT_Dx_resources_"+str(Dx_resources_v)+".log"
output_Dx_INI="GIT_Dx_INI_"+str(Dx_INI_v)+".log"

## Log git commits of used repo's
os.system(" date >> "+str(output_IAP))
os.system(" date >> "+str(output_Chromate))
os.system(" date >> "+str(output_Dx_resources))
os.system(" date >> "+str(output_Dx_INI))
## Print used folder in log files
os.system("echo "+str(IAP_folder)+" >> "+str(output_IAP))
os.system("echo "+str(Chromate_folder)+" >> "+str(output_Chromate))
os.system("echo "+str(Dx_resources_folder)+" >> "+str(output_Dx_resources))
os.system("echo "+str(Dx_INI_folder)+" >> "+str(output_Dx_INI))
## Print GIT version in log files
os.system("echo "+str(IAP_v)+" >> "+str(output_IAP))
os.system("echo "+str(Chromate_v)+" >> "+str(output_Chromate))
os.system("echo "+str(Dx_resources_v)+" >> "+str(output_Dx_resources))
os.system("echo "+str(Dx_INI_v)+" >> "+str(output_Dx_INI))
## Print full GIT log in log files
os.system("git --git-dir="+str(IAP_folder)+".git log >> "+str(output_IAP))
os.system("git --git-dir="+str(Chromate_folder)+".git log >> "+str(output_Chromate))
os.system("git --git-dir="+str(Dx_resources_folder)+".git log >> "+str(output_Dx_resources))
os.system("git --git-dir="+str(Dx_INI_folder)+".git log >> "+str(output_Dx_INI))


# Make run_stats.txt file from all flagstats
os.system("/hpc/cog_bioinf/diagnostiek/production/Dx_resources/get_stats_from_flagstat.pl >run_stats.txt")

sys.exit("Finished")

