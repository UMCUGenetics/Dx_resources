#! /usr/bin/env python
import sys, os, re
import commands
import getopt
from optparse import OptionParser
from optparse import OptionGroup

## Author: M.G. Elferink
## Date: 27-09-2016
## Purpose: Perform re-analysis of the ExonCov tool within a run folder.




if __name__ == "__main__":
        parser = OptionParser();
        group = OptionGroup(parser, "Main options")
	group.add_option("-e", default="/hpc/local/CentOS7/cog_bioinf/ExonCov/ExonCov.py", dest="EXONCALLCOV_PATH", help="EXONCALLCOV_PATH including executable")
	group.add_option("-b", default="/hpc/cog_bioinf/diagnostiek/production/Dx_resources/Tracks/ENSEMBL_UCSC_merged_collapsed_sorted_v2_20bpflank.bed", dest="EXONCALLCOV_BED", help="EXONCALLCOV_BED")
	group.add_option("-g", default="/hpc/cog_bioinf/diagnostiek/production/Dx_resources/Exoncov/gpanels.txt", dest="EXONCALLCOV_PANEL", help="EXONCALLCOV_PANEL")
	group.add_option("-p", default="/hpc/cog_bioinf/diagnostiek/production/Dx_resources/Exoncov/Preferred_transcript_list.txt", dest="EXONCALLCOV_PREF", help="EXONCALLCOV_PREF")
	group.add_option("-n", default="/hpc/cog_bioinf/diagnostiek/production/Dx_resources/Exoncov/NM_ENSEMBL_HGNC.txt", dest="EXONCALLCOV_ENS", help="EXONCALLCOV_ENS")
	group.add_option("-t", default="02:00:00", dest="EXONCALLCOV_TIME", metavar="[STRING]", help="EXONCALLCOV_TIME")
	group.add_option("-q", default="all.q", dest="EXONCALLCOV_QUEUE", metavar="[STRING]", help="EXONCALLCOV_QUEUE")
	group.add_option("-k", default="dedup.realigned.bam$", dest="search_string", metavar="[STRING]", help="search_string [default = dedup.realigned.bam$")
	group.add_option("-m", default="off", dest="EXONCALLCOV_MEM", metavar="[STRING]", help="EXONCALLCOV_MEM")
	group.add_option("-i", default="off", dest="use_ini", action="store_false", help="use information from INI file: OFF [default] or ON")
        group.add_option("-f", default="/hpc/cog_bioinf/diagnostiek/production/Dx_resources/", dest="Dx_resources_folder", metavar="[STRING]", help="Dx_resources_folder")
        group.add_option("-o", default="off", dest="use_custom_folder", help="custom output foldername [default=off]")
	parser.add_option_group(group)
        (opt, args) = parser.parse_args()

	use_ini=opt.use_ini
	search_string= opt.search_string
	Dx_resources_folder=opt.Dx_resources_folder

if use_ini == "ON":
	"Use INI settings"
	try:
		ini=commands.getoutput("find -iname UMCU_DX_*.ini")
		if len(ini.split()) ==1 :
			in_file=open(str(ini.split()[0]),"r").readlines()
		else:
			sys.exit("INI not found or multiple INI found")
	except:
		sys.exit("INI not found")
	
	for line in in_file:
		splitline=line.split()
		if len(splitline)>1:
			if splitline[0]=="EXONCALLCOV_PATH":
				EXONCALLCOV_PATH=splitline[1]
			if splitline[0]=="EXONCALLCOV_BED":
		        	EXONCALLCOV_BED=splitline[1]
			if splitline[0]=="EXONCALLCOV_PREF":
		        	EXONCALLCOV_PREF=splitline[1]
			if splitline[0]=="EXONCALLCOV_PANEL":
		        	EXONCALLCOV_PANEL=splitline[1]
		        if splitline[0]=="EXONCALLCOV_ENS":
		                EXONCALLCOV_ENS=splitline[1]
			if splitline[0]=="EXONCALLCOV_QUEUE":
		                EXONCALLCOV_QUEUE=splitline[1]
			if splitline[0]=="EXONCALLCOV_TIME":
		                EXONCALLCOV_TIME=splitline[1]
			if splitline[0]=="EXONCALLCOV_MEM":
		                EXONCALLCOV_MEM=splitline[1]

elif use_ini== "off":
	print "INI settings not used"
	EXONCALLCOV_PATH=opt.EXONCALLCOV_PATH
	EXONCALLCOV_BED=opt.EXONCALLCOV_BED
	EXONCALLCOV_PREF=opt.EXONCALLCOV_PREF
	EXONCALLCOV_PANEL=opt.EXONCALLCOV_PANEL
	EXONCALLCOV_ENS=opt.EXONCALLCOV_ENS
	EXONCALLCOV_QUEUE=opt.EXONCALLCOV_QUEUE
	EXONCALLCOV_TIME=opt.EXONCALLCOV_TIME
	EXONCALLCOV_MEM=opt.EXONCALLCOV_MEM


EXONCALLCOV_FOLDER="/".join(EXONCALLCOV_PATH.split("/")[0:-1])  ## EXONCALLCOV_PATH is including script file
ExonCov_v=commands.getoutput("git --git-dir="+str(EXONCALLCOV_FOLDER)+"/.git describe --tags")
Dx_resources_v=commands.getoutput("git --git-dir="+str(Dx_resources_folder)+"/.git describe --tags")
if str(opt.use_custom_folder) == "off":
	output_folder="Exoncov_reanalysis_GIT_Dx_resources_"+str(Dx_resources_v)
	output_Dx_resources="Exoncov_reanalysis_GIT_Dx_resources_"+str(Dx_resources_v)+".log"
	output_ExonCov="Exoncov_reanalysis_GIT_ExonCov_"+str(ExonCov_v)+".log"
else:	
	output_folder=str(opt.use_custom_folder)        
	output_Dx_resources=str(output_folder)+"_GIT_Dx_resources_"+str(Dx_resources_v)+".log"
	output_ExonCov=str(output_folder)+"_GIT_ExonCov_"+str(ExonCov_v)+".log"


## Perform ExonCov tool
action="python "+str(EXONCALLCOV_PATH)+ " --queue "+str(EXONCALLCOV_QUEUE)+" -a " +str(EXONCALLCOV_TIME)+ " -c "+ str(EXONCALLCOV_MEM) +" -b "+str(EXONCALLCOV_BED)+" -n "+str(EXONCALLCOV_ENS)+" -p "+ str(EXONCALLCOV_PREF)+ " -l "+str(EXONCALLCOV_PANEL)+ " -d "+ str(search_string) + " -o " +str(output_folder)
print action
os.system(action)

## Write Git-log files
os.system(" date >> "+str(output_ExonCov))
os.system(" date >> "+str(output_Dx_resources))
os.system("echo "+str(ExonCov_v)+" >> "+str(output_ExonCov))
os.system("echo "+str(Dx_resources_v)+" >> "+str(output_Dx_resources))
os.system("echo "+str(ExonCov_v)+" >> "+str(output_ExonCov))
os.system("echo "+str(Dx_resources_v)+" >> "+str(output_Dx_resources))
os.system("git --git-dir="+str(EXONCALLCOV_FOLDER)+"/.git log >> "+str(output_ExonCov))
os.system("git --git-dir="+str(Dx_resources_folder)+"/.git log >> "+str(output_Dx_resources))


if str(opt.use_custom_folder) == "off": 
	action2= "mv Exoncov_reanalysis_GIT*.log "+str(output_folder)
else:
	action2= "mv "+str(output_folder)+"_GIT*.log "+str(output_folder)
action3=" sed -i 's/\./\,/g' "+str(output_folder)+"/*tsv"
action4=" sed -i 's/\./\,/g' "+str(output_folder)+"/*/*tsv"
os.system(action2)
os.system(action3)
os.system(action4)

print "finished"




