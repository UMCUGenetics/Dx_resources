#! /usr/bin/env python
import sys,os,re, commands
# AUTHOR:       M.G. Elferink
# DATE:         26-10-2012
# GC correction for all pickle files

bin_folder=sys.argv[1]
folder= sys.argv[2]
lib = os.listdir(str(folder))
#bin_folder=sys.argv[0].split("/")
#bin_folder= "/".join(bin_folder[0:-1])+"/"

for item in lib:
	if ".out" in item and ".txt" not in item:
		action ="python "+str(bin_folder)+"plot.py "+str(folder)+str(item)+" "+str(folder)+str(item)[0:-4]+ " > "+str(folder)+ str(item)[0:-4]+ "_plot_setting.txt"
		print action
		os.system(action)

sys.exit("Finished")

