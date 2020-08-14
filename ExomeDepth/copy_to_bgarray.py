#! /usr/bin/env python
import sys
import os
import argparse

def copybgarray(args):
    today = date.today().strftime("%d%m%y")
    lines = open(args.run_file,"r").readlines()
    for line in lines:
        splitline = line.split()
        input_folder = "{input_path}/{run}".format(input_path=args.inputfolder, run=splitline[0]) 
        output_folder = "{output_path}/{run}".format(output_path=args.outputfolder, run=splitline[0])

        """Remove exomedepth folder original folder on bgarray"""
        action = "rm -r {output_folder}/exomedepth".format(output_folder=output_folder)
        os.system(action)
        action = "rm -r {output_folder}/ExomeDepth".format(output_folder=output_folder)
        os.system(action)

        """Copy new exomedepth folder into originale folder on bgarray""" 
        action = "rsync -rahu --progress --prune-empty-dirs {input_folder}/exomedepth {output_folder}/".format(input_folder=input_folder, output_folder=output_folder)
        os.system(action)

        """Add new analysis stamp to logbook.txt"""
        action = "echo \"ME {today}\tExomeDepth re-analysis using refset {refset}\" >> {output_folder}/logbook.txt".format(today=today, refset=splitline[1], output_folder=output_folder, run=splitline[0])
        os.system(action)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Input folder IAP/NF analysis')
    parser.add_argument('outputfolder', help='Output folder for new analysis')
    parser.add_argument('run_file', help='Tab delimited file with analysisID + reference set to be used')
    parser.set_defaults(func = copybgarray)
    args = parser.parse_args()
    args.func(args)
