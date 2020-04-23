#! /usr/bin/env python3
import os
import argparse
import subprocess
import settings
from multiprocessing import Pool
import sys

def process(bam):
    bamfile = bam.rstrip("/").split("/")[-1]
    sampleid = bamfile.split("_")[0]
    os.system("mkdir -p {output}/{sample}".format(output = args.outputfolder, sample = sampleid))
    os.system("ln -sd {bam} {output}/{sample}/{bamfile}".format(bam = bam, bamfile = bamfile, sample = sampleid, output = args.outputfolder))
    os.system("ln -sd {bam}.bai {output}/{sample}/{bamfile}.bai".format(bam = bam, bamfile = bamfile,sample = sampleid, output = args.outputfolder))
    os.chdir("{output}/{sample}".format(output = args.outputfolder, sample = sampleid))
    action = ("python {exomedepth} callcnv {output}/{sample} {inputbam} {run} {sample} {refset} --batch".format(
            exomedepth = args.exomedepth,
            output = args.outputfolder,
            inputbam = bamfile,
            run = args.inputfolder.rstrip("/").split("/")[-1],
            sample = bam.split("/")[-1].split("_")[0],
            refset = args.refset
        ))
    os.system(action)

    os.chdir("{output}".format(output = args.outputfolder))
    os.system("mkdir -p {output}/UMCU/{sampleid}".format(output = args.outputfolder, sampleid = sampleid))
    os.system("mkdir -p {output}/HC/{sampleid}".format(output = args.outputfolder, sampleid = sampleid))
    os.system("mkdir -p {output}/logs".format(output = args.outputfolder))
    os.system("mkdir -p {output}/igv_tracks".format(output = args.outputfolder))
    os.system("mv {output}/{sampleid}/*log {output}/logs/".format(sampleid = sampleid, output = args.outputfolder))
    os.system("mv {output}/{sampleid}/*igv {output}/igv_tracks/".format(sampleid = sampleid, output = args.outputfolder))
    os.system("mv {output}/{sampleid}/HC*vcf {output}/HC/".format(sampleid = sampleid, output = args.outputfolder))
    os.system("mv {output}/{sampleid}/UMCU*vcf {output}/UMCU/".format(sampleid = sampleid, output = args.outputfolder))
    os.system("mv {output}/{sampleid}/HC* {output}/HC/{sampleid}".format(sampleid = sampleid, output = args.outputfolder))
    os.system("mv {output}/{sampleid}/UMCU* {output}/UMCU/{sampleid}".format(sampleid = sampleid, output = args.outputfolder))
    os.system("rm -r {output}/{sample}".format(output = args.outputfolder, sample = sampleid))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Path to input folder containing BAM files to process')
    parser.add_argument('outputfolder', help='Path to output folder')
    parser.add_argument('simjobs', help='number of simultanious samples to proces. Note: make sure similar threads are reseved in session!')
    parser.add_argument('--refset', default = settings.refset, help='Reference set to be used')
    parser.add_argument('--exomedepth', default = "/hpc/diaggen/software/production/Dx_resources/ExomeDepth/run_ExomeDepth.py", help='Full path to exomedepth script')
    args = parser.parse_args()

    bams = subprocess.getoutput("find -L {0} -iname \"*.realigned.bam\"".format(args.inputfolder)).split()
    print("Number of BAM files = "+str(len(bams)))
   
    with Pool(processes=int(args.simjobs)) as pool:
        result = pool.map(process, bams, 1)
