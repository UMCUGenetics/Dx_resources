#! /usr/bin/env python3
import os
import argparse
import subprocess
import glob
from multiprocessing import Pool
import settings
from datetime import date

import sys

def process(bam):
    bamfile = bam.rstrip("/").split("/")[-1]
    if args.pipeline == "iap":
        sampleid = bamfile.split("_")[0]
    elif args.pipeline == "nf":
        sampleid = bamfile.split(".")[0]

    if args.refsetlist:
        if sampleid in refset_dic:
            refset = refset_dic[sampleid]
        else:
            refset = args.refset
    else:
        refset = args.refset

    os.system("mkdir -p {output}/{sample}".format(output = args.outputfolder, sample = sampleid))
    os.system("ln -sd {bam} {output}/{sample}/{bamfile}".format(bam = bam, bamfile = bamfile, sample = sampleid, output = args.outputfolder))
    os.system("ln -sd {bam}.bai {output}/{sample}/{bamfile}.bai".format(bam = bam, bamfile = bamfile, sample = sampleid, output = args.outputfolder))
    os.chdir("{output}/{sample}".format(output = args.outputfolder, sample = sampleid))

    if args.pipeline == "iap":
        action = "python {exomedepth} callcnv {output}/{sample} {inputbam} {run} {sample} {refset} --expectedCNVlength {length} --pipeline iap".format(
            exomedepth = args.exomedepth,
            output = args.outputfolder,
            inputbam = bamfile,
            run = args.inputfolder.rstrip("/").split("/")[-1],
            sample = sampleid,
            refset = refset,
            length = args.expectedCNVlength
        )
    elif args.pipeline == "nf":
        action = "python {exomedepth} callcnv {output}/{sample} {inputbam} {run} {sample} {refset} --expectedCNVlength {length} --pipeline nf".format(
            exomedepth = args.exomedepth,
            output = args.outputfolder,
            inputbam = bamfile,
            run = args.inputfolder.rstrip("/").split("/")[-1],
            sample = sampleid,
            refset = refset,
            length = args.expectedCNVlength
        )

    os.system(action)

    os.chdir("{output}".format(output = args.outputfolder))
    os.system("mkdir -p {output}/logs".format(output = args.outputfolder))
    os.system("mkdir -p {output}/igv_tracks".format(output = args.outputfolder))
    os.system("mkdir -p {output}/UMCU/".format(output = args.outputfolder))
    os.system("mkdir -p {output}/HC/".format(output = args.outputfolder))
    os.system("mv {output}/{sampleid}/*.xml {output}/".format(sampleid = sampleid, output = args.outputfolder)) 
    os.system("mv {output}/{sampleid}/*.log {output}/logs/".format(sampleid = sampleid, output = args.outputfolder))
    os.system("mv {output}/{sampleid}/*.igv {output}/igv_tracks/".format(sampleid = sampleid, output = args.outputfolder))
    os.system("mv {output}/{sampleid}/HC*.vcf {output}/HC/".format(sampleid = sampleid, output = args.outputfolder))
    os.system("mv {output}/{sampleid}/UMCU*.vcf {output}/UMCU/".format(sampleid = sampleid, output = args.outputfolder))
    os.system("rm -r {output}/{sample}".format(output = args.outputfolder, sample = sampleid))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Path to root folder of analysis')
    parser.add_argument('outputfolder', help='Path to output folder')
    parser.add_argument('simjobs', help='number of simultaneous samples to proces. Note: make sure similar threads are reseved in session!')
    parser.add_argument('--pipeline', default='nf', choices=['nf', 'iap'], help='pipeline used for sample processing (nf = nexflow (default), IAP = illumina analysis pipeline)')
    parser.add_argument('--refset', default = settings.refset, help='Reference set to be used')
    parser.add_argument('--refsetlist', help='Tab delimited file with SampleID and RefsetID to be used. If samples are not present in the file, default refset is used')
    parser.add_argument('--exomedepth', default = "/hpc/diaggen/software/production/Dx_resources/ExomeDepth/run_ExomeDepth.py", help='Full path to exomedepth script')
    parser.add_argument('--expectedCNVlength',default=settings.expectedCNVlength, help='expected CNV length (basepairs) taken into account by ExomeDepth [default expectedCNVlength in settings.py]')
    args = parser.parse_args()

    today = date.today().strftime("%d%m%y")
    user = subprocess.getoutput("whoami") 

    """ Check if previous exomedepth analysis is already present """
    if os.path.isdir(args.outputfolder):
        print("Run folder not empty: assuming exomedepth has been runned. Current exomedepth folder will be archived")
        """Check is archive folder exists. If this is the case, relative path in IGV session should not be changed."""
        archivefolder = False
        if os.path.isdir("{outputfolder}/archive_{today}/exomedepth".format(outputfolder=args.outputfolder, today=today)):
            archivefolder = True
        """ Make archive folder"""
        os.system("mkdir -p {outputfolder}/archive_{today}/exomedepth".format(outputfolder=args.outputfolder, today=today))
        """ Move original data to archive folder """
        os.system("mv {outputfolder}/* {outputfolder}/archive_{today}/exomedepth/".format(outputfolder=args.outputfolder, today=today))
        """ Rename relative paths in IGV sessions"""
        if archivefolder == False:
            os.system("sed -i 's/\"\.\.\//\"\.\.\/\.\.\/\.\.\//g' {outputfolder}/archive_{today}/*/*xml".format(outputfolder=args.outputfolder, today=today))
        """ Check if archive folder were already present in archived folder, and move the to the correct location."""
        if glob.glob("{outputfolder}/archive_{today}/exomedepth/archive*".format(outputfolder=args.outputfolder, today=today)):
            os.system("mv {outputfolder}/archive_{today}/exomedepth/archive* {outputfolder}/".format(outputfolder=args.outputfolder, today=today))       

    """Find BAM files to be processed"""
    if args.pipeline == "iap":
        bams = set(glob.glob("{}/**/*.realigned.bam".format(args.inputfolder), recursive=True)) - set(glob.glob("{}/[eE]xome[dD]epth*/**/*.realigned.bam".format(args.inputfolder), recursive=True))
    elif args.pipeline == "nf":
        bams = glob.glob("{}/bam_files/**/*.bam".format(args.inputfolder), recursive=True)  
    print("Number of BAM files = "+str(len(bams)))

    if len(bams) > 0: 
        os.system("echo \"{user} {today}\tExomeDepth reanalysis performed\" >> {inputfolder}/logbook.txt".format(user=user, today=today, inputfolder=args.inputfolder))
    else:
        sys.exit("no bam files detected")

    if args.refsetlist:
        refset_dic = {}
        refset_lines = open (args.refsetlist,"r").readlines()
        for line in refset_lines:
           splitline = line.split()
           if splitline[0] not in refset_dic:
               refset_dic[splitline[0]] = splitline[1]

    """Start exomedepth re-analysis"""
    with Pool(processes=int(args.simjobs)) as pool:
        result = pool.map(process, bams, 1)
