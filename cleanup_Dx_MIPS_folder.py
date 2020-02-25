#!/usr/bin/env python
import commands
import fnmatch
import glob
import os
import shutil


# Create git logfiles for files in inputlist
def logging(inputlist, folderlist):
    for repo in inputlist:
        index = inputlist.index(repo)
        repo_folder = folderlist[index]
        repo_v = commands.getoutput("git --git-dir={}.git describe --tags".format(repo_folder))
        output_repo = "GIT_{}_{}.log".format(repo, repo_v)
        os.system(" date >> "+str(output_repo))
        os.system("echo "+str(repo_folder)+" >> "+str(output_repo))
        os.system("echo "+str(repo_v)+" >> "+str(output_repo))
        os.system("git --git-dir="+str(repo_folder)+".git log >> "+str(output_repo))


# Move disapproved samples to bad data folder
def move(vcf):
    for i in os.listdir('.'):
        if i.startswith(vcf):
            shutil.move(i, badfolder)


# Variables
inputlist = ["IAP", "Dx_resources", "Dx_INI", "Dx_tracks", "mips", "fingerprinting"]
folderlist = [
    "/hpc/diaggen/software/production/IAP",
    "/hpc/diaggen/software/production/Dx_resources/",
    "/hpc/diaggen/software/production/Dx_INI/",
    "/hpc/diaggen/software/production/Dx_tracks/",
    "/hpc/diaggen/software/production/mips/",
    "/hpc/diaggen/software/production/fingerprinting/"
]

badfolder = 'disapprovedVCFs'

# Create GIT log files
print '\nMaking GIT log files...'

logging(inputlist, folderlist)

# Check results

print "\nCheck results..."

os.chdir('./fingerprint')
os.mkdir(badfolder)
for item in os.listdir('.'):
    if fnmatch.fnmatch(item, '*.vcf'):
        with open(item, 'r') as f:
            refcov = 0  # Total reference reads for all homozygous (1/1) calls
            totalcov = 0  # Total coverage for all homozygous calls
            homaltcount = 0  # Number of homozygous calls
            ycount = 0  # Sum of coverage for two Y SNPs
            lowcovcount = 0  # Number of SNPs with <15X coverage
            disbalancecount = 0  # Number of heterozygous (0/1) calls with allelefrequency <0.2 or >0.8

            for line in f:
                if not line.startswith('#'):
                    line = line.split()

                    # Parse Genotype format
                    gt_format = line[8].split(':')
                    gt_index = gt_format.index('GT')

                    # Parse sample genotype
                    gt_values = line[9].split(':')
                    #gt_values = gt_values.replace(':', '\t').replace(',', '\t')
                    gt_value = gt_values[gt_index]

                    if line[0] == 'Y':
                        if gt_value != './.':
                            ycount += int(gt_values[gt_format.index('DP')])
                    elif gt_value == '1/1':
                        homaltcount += 1
                        if int(gt_values[gt_format.index('DP')]) < 15:
                            lowcovcount += 1
                        refcov += int(gt_values[gt_format.index('AD')].split(',')[0])
                        totalcov += int(gt_values[gt_format.index('DP')])
                    elif gt_value != '1/1':
                        if gt_value == './.':
                            lowcovcount += 1
                        elif gt_value == '0/0':
                            if int(gt_values[gt_format.index('DP')]) < 15:
                                lowcovcount += 1
                        else:
                            if int(gt_values[gt_format.index('DP')]) < 15:
                                lowcovcount += 1
                    if gt_value == '0/1':
                        af_value = int(gt_values[gt_format.index('AD')].split(',')[0]) / float(int(gt_values[gt_format.index('DP')]))
                        if af_value > 0.8 or af_value < 0.2:
                            disbalancecount += 1

            contamination = refcov / float(totalcov)

            result = item, str(lowcovcount), str(homaltcount), str(round(contamination, 6)), item[8], str(ycount), str(disbalancecount)
            resultprint = '\t'.join(result)
            os.system("echo " + resultprint + " >> ../logs.tmp")

            with open('../message.tmp', 'a') as message:
                gender = ['M', 'F', 'O']
                if result[4] not in gender:
                    message.write('### {}: report filename issue to lab'.format(item) + '\n')
                if int(result[1]) > 15:
                    message.write('### {}: >15 SNPs with <15X coverage ({}) --> disapproved'.format(item, result[1]) + '\n')
                    move(item)
                elif int(result[6]) > 8:
                    message.write('### {}: >8 heterozygous SNPs with <20% MAF ({}) --> disapproved'.format(item, result[6]) + '\n')
                    move(item)
                elif int(result[2]) < 8:
                    message.write('### {}: <8 homozygous ALT SNPs called ({}) --> disapproved'.format(item, result[2]) + '\n')
                    move(item)
                elif result[4] == 'F' and int(result[5]) > 100 or result[4] == 'M' and int(result[5]) < 100:
                    message.write('### {}: gender {} with {} reads on chromosome Y, discuss with lab and disapprove if needed'.format(item, result[4], result[5]) + '\n')

# Write logbook
print('\nCreating logbook...')

os.chdir('..')
os.system("cat logs.tmp | sort | sed 's/ /\t/g' > logbook.tmp")
os.system("cat message.tmp | sort > messagelog.tmp")
with open('logbook.txt', 'a') as out:
    with open('messagelog.tmp', 'r') as message, open('logbook.tmp', 'r') as logs:
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
