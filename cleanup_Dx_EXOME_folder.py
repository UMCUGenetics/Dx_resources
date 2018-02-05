#! /usr/bin/env python
import sys, os, re
import commands
import getopt
import fnmatch
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', action='store_true', default=False, dest='overrule', help='clean-up data even if kinship is not in agreement with pedigree. Default= off')
args = parser.parse_args()

kinship_setting=[0.177,0.354]

def pedfile_dic(pwd):
    """Create ped file dictonary."""
    for file in os.listdir(pwd):
        if fnmatch.fnmatch(file, '*[Pp][eE][dD]'):
	    ped_file = '{}/{}'.format(pwd,file)
    lines = open(ped_file, "r").readlines()
    relation = {}
    family = {}
    for line in lines:
        splitline = line.split()
        try:
            family[splitline[1]]
	except:
	    family[splitline[1]] = splitline[0]

        if splitline[2] == "0" and splitline[3] == "0":	 # both parents are unknown
            pass
        else:
            if str(splitline[2]) != "0": 
                relation[str(splitline[1]) + "_" + str(splitline[2])] = "parent"
                relation[str(splitline[2]) + "_" + str(splitline[1])] = "parent"
            if str(splitline[3]) != "0":
                relation[str(splitline[1]) + "_" + str(splitline[3])] = "parent"
                relation[str(splitline[3]) + "_" + str(splitline[1])] = "parent"
	    if str(splitline[3]) != "0":
	        relation[str(splitline[3])] = "mother"
	    if str(splitline[2]) != "0":	
                relation[str(splitline[2])] = "father"
    return relation, family

def check_kinship(kinship, pedigree):
    """Check kinship based on pedigree information."""
    kin_warn = []
    if "ID1" not in kinship[0]:
        sample1 = kinship[0]
        sample2 = kinship[1]
        relation = pedigree[0]
        family = pedigree[1]
        try:  # sample relation is parent-child
            relation[sample1+"_"+sample2]
            if family[sample1] == family[sample2] and relation[sample1+"_"+sample2] == "parent":
                if float(kinship[2]) > float(kinship_setting[0]) and float(kinship[2]) < float(kinship_setting[1]):
                    pass
                else:
                    kin_warn += ["PARENT-CHILD relation is not correct "+" ".join(kinship)]
        except:  # sample relation is not parent-child
            if family[sample1] == family[sample2]:  # either parents or sibling
                try:
                     # parent-parent
                     relation[sample1]
                     relation[sample2]
                     #if float(kinship[2]) > float(kinship_setting[0]) and float(kinship[2]) < float(kinship_setting[1]):
                     if float(kinship[2]) > float(kinship_setting[0]): # if parent-parent is more than lower-threshold
                         kin_warn += ["PARENT-PARENT relation is not correct (parent-child or sibling relationship) "+" ".join(kinship)]
                except:
                    # siblings
                    if float(kinship[2]) > float(kinship_setting[0]) and float(kinship[2]) < float(kinship_setting[1]):
                        pass  # kinship is correct
                    else:
                        kin_warn += ["SIBLING-SIBLING relation is not correct "+" ".join(kinship)]
            else:   # unrelated samples
                #if float(kinship[2]) > float(kinship_setting[0]) and float(kinship[2]) < float(kinship_setting[1]):
                if float(kinship[2]) > float(kinship_setting[0]): # if unrelated is more than lower-threshold
                    kin_warn += ["UNRELATED relation is not correct (kinship is >0.177) "+" ".join(kinship)]
    return kin_warn

def extract_qc(pwd):
    kinship_file = '{}/{}.kinship'.format(pwd, pwd.split("/")[-1])
    kinship_lines = open(str(kinship_file), "r").readlines()
    warning=[]
    pedigree = pedfile_dic(pwd)
    for line in kinship_lines:
        splitline = line.split()
        printline = [splitline[1], splitline[3], splitline[7]]
        warning += check_kinship(printline, pedigree)
    return warning

print "Calculating run stats from flagstat files\n" 
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

# Check kinship
pwd= commands.getoutput("pwd")

warning= extract_qc(pwd)
if len(warning) > 0:
    print "#"*100, "\n","#"*100, "\n"," "*30,"WARNING!!\n","#"*100,"\n","#"*100, "\n"
    for line in warning:
        print line
    if args.overrule == True:
        print "\n"+"#"*100,"\nWarning over ruled with -f command"
    else:
        sys.exit("\n"+"#"*100+"\nclean-up not completed!"+"\n"+"#"*100)
else:
    print "Kinship checked with pedigree: ok!\n"

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

print "Making GIT log files\n"
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

print "Moving files\n" 
# move unused VCFs in Redundant folder
os.system("mkdir Redundant_VCF_files/")
#os.system("find -iname \"*phased.vcf*\" -exec mv {} Redundant_VCF_files/ \;")
os.system("find -iname \"*raw_variants.vcf*\" -exec mv {} Redundant_VCF_files/ \;")

print "Creating .bam.bai files\n"
# cp .bai into bam.bai
pwd= commands.getoutput("pwd")
action="find "+str(pwd)+" -iname \"*dedup.realigned.bai\"" 
files=commands.getoutput(action).split()
for file in files:
	action= "cp " + file+ " "+ file[:-1]+"m.bai"
	print action
	os.system(action)

print "Converting Exoncov files\n"
# change dots in Exoncov files into comma's
action2=" sed -i 's/\./\,/g' "+str(pwd)+"/Exoncov*/*tsv"	## include use Exoncov output folder?
action3=" sed -i 's/\./\,/g' "+str(pwd)+"/Exoncov*/*/*tsv"	## include use Exoncov output folder?
os.system(action2)
os.system(action3)

print "Finished"

