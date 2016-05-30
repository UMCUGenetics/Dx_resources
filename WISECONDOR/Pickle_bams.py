#! /usr/bin/env python
import sys,os,re, commands
# AUTHOR:       M.G. Elferink
# DATE:         26-10-2012
# Make pickle files for all BAM files within folder

bin_folder=sys.argv[1]
folder= sys.argv[2]
#bin_folder=sys.argv[0].split("/")
#bin_folder= "/".join(bin_folder[0:-1])+"/"
lib = os.listdir(str(folder))

for item in lib:
	if ".bam" in item and "bai" not in item:
		action= "/hpc/local/CentOS7/cog_bioinf/sambamba_v0.5.8/sambamba view "+str(folder)+str(item) +" -t 4 --filter=\"mapping_quality >= 40 and [NM]<=1 and [SA] == null and not duplicate and not failed_quality_control \"| python "+str(bin_folder)+"consam.py -outfile "+str(folder) + str(item)[0:-4]+ ".pickle"
		print action
		os.system(action)

sys.exit("Finished")

