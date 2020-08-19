#! /usr/bin/env python
import os
from datetime import date
import os.path
import argparse

def reanalyse(args):
    today = date.today().strftime("%d%m%y")
    lines = open(args.run_file,"r").readlines()
    for line in lines:
        splitline = line.split()
        input_folder = "{input_path}/{run}".format(input_path=args.inputfolder, run=splitline[0])
        output_folder = "{output_path}/{run}/exomedepth/".format(output_path=args.outputfolder, run=splitline[0])

        os.system("mkdir -p {output_folder}/logs/".format(output_folder=output_folder))
        os.system("mkdir -p {output_folder}/archive_{today}/".format(output_folder=output_folder, today=today))

        """Copy original exomedepth data into archive folder"""
        if os.path.isdir("{input_folder}/ExomeDepth".format(input_folder=input_folder)) or os.path.isdir("{input_folder}/exomedepth".format(input_folder=input_folder)):
            os.system("rsync -rahuL --progress {input_folder}/[Ee]xome[Dd]epth {output_folder}/archive_{today}/".format(input_folder=input_folder, output_folder=output_folder, today=today))

        """Adjust relatieve paths in XML files (2 extra ../): only for the file that are the latest in the current storage"""
        if os.path.isdir("{output_folder}/archive_{today}/ExomeDepth".format(output_folder=output_folder, today=today)) or os.path.isdir("{output_folder}/archive_{today}/exomedepth".format(output_folder=output_folder, today=today)):
            os.system("sed -i 's/\"\.\.\//\"\.\.\/\.\.\/\.\.\//g' {output_folder}/archive_{today}/*/*xml".format(output_folder=output_folder, today=today))
        
        """ Check if archive folder within new archive folder and mv this to root folder of new exomedepth folder. 
        This will prevent recursive archive folders, and prevent problems with relative paths."""
        if glob.glob("{output_folder}/archive_{today}/exomedepth/archive*".format(output_folder=output_folder, today=today)): 
            os.system("mv {output_folder}/archive_{today}/exomedepth/archive* {output_folder}/".format(output_folder=output_folder, today=today))

        """Perform new ExomeDepth analysis """
	os.chdir("{output_folder}/".format(output_folder=output_folder))
        action = "sbatch -t 6:0:00 --mem=80G -c 16 --mail-type=FAIL --mail-user=m.elferink@umcutrecht.nl --wrap=\"source /hpc/diaggen/software/production/Dx_resources/ExomeDepth/venv/bin/activate && /hpc/diaggen/software/production/Dx_resources/ExomeDepth/submit_batch_exomedepth.py {input_folder}/ {output_folder} 7 --pipeline {pipeline} --refset {refset}\"".format(input_folder=input_folder, output_folder=output_folder, pipeline=args.pipeline, refset=splitline[1])
        os.system(action)

        """Copy logbook.txt to new folder, and add new analysis"""
        if os.path.isfile("{input_folder}/logbook.txt".format(input_folder=input_folder)):
            os.system("rsync -rahuL --progress {input_folder}/logbook.txt {output_folder}/{run}/".format(input_folder=input_folder, output_folder=args.outputfolder, run=splitline[0]))
        os.system("echo \"{user} {today}\tExomeDepth reanalysis using refset {refset}\" >> {output_path}/{run}/logbook.txt".format(user=args.user, today=today, refset=splitline[1], output_path=args.outputfolder, run=splitline[0]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Input folder IAP/NF analysis')
    parser.add_argument('outputfolder', help='Output folder for new analysis')
    parser.add_argument('run_file', help='Tab delimited file with analysisID + reference set to be used')
    parser.add_argument('pipeline', choices=['nf', 'iap'], help='pipeline used for sample processing (nf = nexflow, iap = illumina analysis pipeline')
    parser.add_argument('user', help='user abbreviation. Will be mentioned in logbook.txt')
    parser.set_defaults(func = reanalyse)
    args = parser.parse_args()
    args.func(args)



