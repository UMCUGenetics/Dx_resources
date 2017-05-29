#! /usr/bin/env python
import os
import commands

# Input folders
IAP_folder="/hpc/local/CentOS7/cog_bioinf/IAP_Dx/"
Dx_resources_folder="/hpc/cog_bioinf/diagnostiek/production/Dx_resources/"
Dx_INI_folder="/hpc/cog_bioinf/diagnostiek/production/Dx_INI/"
Dx_tracks_folder="/hpc/cog_bioinf/diagnostiek/production/Dx_tracks/"
Dx_mips_folder="/hpc/cog_bioinf/diagnostiek/production/Dx_mips/"
Dx_fingerprinting_folder="/hpc/cog_bioinf/diagnostiek/production/Dx_fingerprinting/"
# Get GIT tag version
IAP_v=commands.getoutput("git --git-dir="+str(IAP_folder)+".git describe --tags")
Dx_resources_v=commands.getoutput("git --git-dir="+str(Dx_resources_folder)+".git describe --tags")
Dx_INI_v=commands.getoutput("git --git-dir="+str(Dx_INI_folder)+".git describe --tags")
Dx_tracks_v=commands.getoutput("git --git-dir="+str(Dx_tracks_folder)+".git describe --tags")
Dx_mips_v=commands.getoutput("git --git-dir="+str(Dx_mips_folder)+".git describe --tags")
Dx_fingerprinting_v=commands.getoutput("git --git-dir="+str(Dx_fingerprinting_folder)+".git describe --tags")
# Output files
output_IAP="GIT_IAP_"+str(IAP_v)+".log"
output_Dx_resources="GIT_Dx_resources_"+str(Dx_resources_v)+".log"
output_Dx_INI="GIT_Dx_INI_"+str(Dx_INI_v)+".log"
output_Dx_tracks="GIT_Dx_tracks_"+str(Dx_tracks_v)+".log"
output_Dx_mips="GIT_Dx_mips_"+str(Dx_mips_v)+".log"
output_Dx_fingerprinting="GIT_Dx_fingerprinting_"+str(Dx_fingerprinting_v)+".log"

print "Making GIT log files\n"

## Log git commits of used repo's
os.system(" date >> "+str(output_IAP))
os.system(" date >> "+str(output_Dx_resources))
os.system(" date >> "+str(output_Dx_INI))
os.system(" date >> "+str(output_Dx_tracks))
os.system(" date >> "+str(output_Dx_mips))
os.system(" date >> "+str(output_Dx_fingerprinting))
## Print used folder in log files
os.system("echo "+str(IAP_folder)+" >> "+str(output_IAP))
os.system("echo "+str(Dx_resources_folder)+" >> "+str(output_Dx_resources))
os.system("echo "+str(Dx_INI_folder)+" >> "+str(output_Dx_INI))
os.system("echo "+str(Dx_tracks_folder)+" >> "+str(output_Dx_tracks))
os.system("echo "+str(Dx_mips_folder)+" >> "+str(output_Dx_mips))
os.system("echo "+str(Dx_fingerprinting_folder)+" >> "+str(output_Dx_fingerprinting))
## Print GIT version in log files
os.system("echo "+str(IAP_v)+" >> "+str(output_IAP))
os.system("echo "+str(Dx_resources_v)+" >> "+str(output_Dx_resources))
os.system("echo "+str(Dx_INI_v)+" >> "+str(output_Dx_INI))
os.system("echo "+str(Dx_tracks_v)+" >> "+str(output_Dx_tracks))
os.system("echo "+str(Dx_mips_v)+" >> "+str(output_Dx_mips))
os.system("echo "+str(Dx_fingerprinting_v)+" >> "+str(output_Dx_fingerprinting))
## Print full GIT log in log files
os.system("git --git-dir="+str(IAP_folder)+".git log >> "+str(output_IAP))
os.system("git --git-dir="+str(Dx_resources_folder)+".git log >> "+str(output_Dx_resources))
os.system("git --git-dir="+str(Dx_INI_folder)+".git log >> "+str(output_Dx_INI))
os.system("git --git-dir="+str(Dx_tracks_folder)+".git log >> "+str(output_Dx_tracks))
os.system("git --git-dir="+str(Dx_mips_folder)+".git log >> "+str(output_Dx_mips))
os.system("git --git-dir="+str(Dx_fingerprinting_folder)+".git log >> "+str(output_Dx_fingerprinting))

print "Finished"
