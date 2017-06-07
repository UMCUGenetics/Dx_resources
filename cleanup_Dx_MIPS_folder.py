#!/usr/bin/env python

import fnmatch
import os
import shutil
import commands

# Versioning input folders
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
# Versioning output files
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

print "Check results:"

# Result folder and files
badfolder = 'disapprovedVCFs'
preoutput = '../tmp.txt'
output = '../logbook.txt'

## Generate results
os.chdir('./fingerprint')
os.mkdir(badfolder)
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*.vcf'):
        with open(file, 'r') as f:
            column2 = 0
            column4 = 0
            hom_count = 0
            y_reads = 0 
            low_cov = 0
            for line in f:
                if not line.startswith('#'):
                    line = line.split()
                    values = line[9]
                    values = values.replace(':', '\t').replace(',', '\t')
                    values = values.split('\t')
                    if line[0] == 'Y':
                        if line[9] != './.':
                            y_reads += int(values[1])
                    elif values[0] == '1/1':
                        if int(values[3]) < 15:
                            low_cov += 1
                        for hom in values[1]:
                            hom_count += 1
                        column2 += int(values[1])
                        column4 += int(values[3])
                    elif values[0] != '1/1':
                        if values[0] == './.':
                            low_cov += 1
                        elif values[0] == '0/0':
                            if int(values[1]) < 15:
                                low_cov += 1
                        else:
                            if int(values[3]) < 15:
                                low_cov += 1
            contamination = column2 / float(column4)
            result = file, str(low_cov), str(hom_count), str(round(contamination,6)), file[8], str(y_reads)            

### Print warnings and replace files
            if result[4] != 'M' and result[4] != 'F':
                print('### For ' + file + ', report filename issue (no M/F) to lab')
            if int(result[1]) > 4: 
                print('### ' + file + ' has ' + str(low_cov) + ' SNPs with coverage <15X --> disapproved!')
                for sample in os.listdir('.'):
                    if sample.startswith(file):
                        shutil.move(sample, badfolder)                
	    elif float(result[3]) > 0.02:
                print('### ' + file + ' has a contamination value of ' + str("{:.2%}".format(contamination)) + ', disapproved!')
                for sample in os.listdir('.'):
                    if sample.startswith(file):
                        shutil.move(sample, badfolder)
            elif int(result[2]) < 10:
                print('### ' + file + ' has ' + str(hom_count) + ' SNPs called homozygous, check VCF for contamination!')
            elif result[4] == 'M' and int(result[5]) < 30:
                print('### ' + file + ' has ' + str(y_reads) + ' reads on chromosome Y, if no low coverage message, check VCF!')
            elif result[4] == 'F' and int(result[5]) > 30:
                print('### ' + file + ' has ' + str(y_reads) + ' reads on chromosome Y, if no contamination message, check VCF!')

#### Create logbook
            str_result = '\t'.join(list(result))
            os.system("echo " + str_result + " >> " + preoutput)

##### Clean-up logbook
f = open(output,'w')
f.write('\n')
f.write('\n')
f.close()
os.system("cat " + preoutput + " | sort | sed 's/ /\t/g' >> " + output)
os.remove(preoutput)

print('\nFinished!')
