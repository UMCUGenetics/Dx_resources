#!/usr/bin/env python

import commands
import fnmatch
import glob
import os
import shutil

# Functions

""" Create git logfiles for files in inputlist """ 
def logging(inputlist):
    for repo in inputlist:
        index = inputlist.index(repo)
        repo_folder = folderlist[index]
        repo_v = commands.getoutput("git --git-dir={}.git describe --tags".format(repo_folder))
        output_repo = "GIT_{}_{}.log".format(repo, repo_v)
        os.system(" date >> "+str(output_repo))
        os.system("echo "+str(repo_folder)+" >> "+str(output_repo))
        os.system("echo "+str(repo_v)+" >> "+str(output_repo))
        os.system("git --git-dir="+str(repo_folder)+".git log >> "+str(output_repo))

""" Move disapproved samples to bad data folder """
def move(vcf):
    for i in os.listdir('.'):
        if i.startswith(vcf):
            shutil.move(i, badfolder)

# Variables

inputlist = ["IAP", "Dx_resources", "Dx_INI", "Dx_tracks", "Dx_mips", "Dx_fingerprinting"]
folderlist = ["/hpc/local/CentOS7/cog_bioinf/IAP_Dx/", 
              "/hpc/cog_bioinf/diagnostiek/production/Dx_resources/",
              "/hpc/cog_bioinf/diagnostiek/production/Dx_INI/",
              "/hpc/cog_bioinf/diagnostiek/production/Dx_tracks/",
              "/hpc/cog_bioinf/diagnostiek/production/Dx_mips/",
              "/hpc/cog_bioinf/diagnostiek/production/Dx_fingerprinting/"]

badfolder = 'disapprovedVCFs'

## Create GIT log files

print "\nMaking GIT log files..."

logging(inputlist)

## Check results

print "\nCheck results..."

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
            false_het = 0
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
                    if values[0] == '0/1':
                        af_value = int(values[1]) / float(int(values[3]))
                        if af_value > 0.8 or af_value < 0.2:
                            false_het += 1             

            contamination = column2 / float(column4)

            result = file, str(low_cov), str(hom_count), str(round(contamination,6)), file[8], str(y_reads)            
            str_result = '\t'.join(list(result))
            os.system("echo " + str_result + " >> ../logs.tmp")

            with open('../message.tmp','a') as message:
                gender = ['M', 'F', 'O']
                if result[4] not in gender:
                    message.write('### {}: REPORT filename issue TO LAB'.format(file) + '\n')
                if int(result[1]) > 4: 
                    message.write('### {}: {} SNPs with <15X coverage --> disapproved'.format(file, result[1]) + '\n') 
                    move(file)                
                elif float(result[3]) > 0.02:
                    message.write('### {}: Contamination value {:.2%} --> disapproved'.format(file, contamination) + '\n')
                    move(file)
                elif int(result[2]) < 10:
                    if false_het > 3: 
                        message.write('### {}: Only {} homozygous SNPs called, {} unbalanced heterozygous calls --> disapproved!'.format(file, hom_count, false_het) + '\n')
                        move(file)
                    else:
                        message.write('### {}: Only {} homozygous SNPs called but no contamination detected --> OK'.format(file, hom_count) + '\n')
                elif result[4] == 'F' and int(result[5]) > 30:
                    if false_het > 3: 
                        message.write('### {}: {} reads on chromosome Y, {} unbalanced heterozygous calls --> disapproved!'.format(file, y_reads, false_het) + '\n')
                        move(file)
                    else:
                        message.write('### {}: {} reads on chromosome Y but no contamination detected, REPORT TO LAB'.format(file, y_reads) + '\n')
                if result[4] == 'M' and int(result[5]) < 100:
                    message.write('### {}: Only {} reads on chromosome Y, REPORT TO LAB'.format(file, y_reads) + '\n')
                
## Write logbook

print('\nCreating logbook...')

os.chdir('..')
os.system("cat logs.tmp | sort | sed 's/ /\t/g' > logbook.tmp")
os.system("cat message.tmp | sort > messagelog.tmp")
with open('logbook.txt','a') as out:
    with open('messagelog.tmp','r') as message, open('logbook.tmp','r') as logs:
        out.write('\n')
        for mline in message:
            mline = mline.rstrip()
            out.write(mline + '\n')
        out.write('\n')
        for lline in logs:
            lline = lline.rstrip()
            out.write(lline + '\n')

for tmp in glob.glob("*.tmp"):
    os.remove(tmp)

print('\nFinished!\n')
