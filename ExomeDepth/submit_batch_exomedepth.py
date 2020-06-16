#! /usr/bin/env python3
import os
import argparse
import subprocess
from multiprocessing import Pool
import settings

def process(bam):
    bamfile = bam.rstrip("/").split("/")[-1]
    if args.pipeline == "iap":
        sampleid = bamfile.split("_")[0]
    elif args.pipeline == "nf":
        sampleid = bamfile.split(".")[0]
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
            refset = args.refset,
            length = args.expectedCNVlength
        )
    elif args.pipeline == "nf":
        action = "python {exomedepth} callcnv {output}/{sample} {inputbam} {run} {sample} {refset} --expectedCNVlength {length} --pipeline nf".format(
            exomedepth = args.exomedepth,
            output = args.outputfolder,
            inputbam = bamfile,
            run = args.inputfolder.rstrip("/").split("/")[-1],
            sample = sampleid,
            refset = args.refset,
            length = args.expectedCNVlength
        )

    os.system(action)

    os.chdir("{output}".format(output = args.outputfolder))
    os.system("mkdir -p {output}/UMCU/{sampleid}".format(output = args.outputfolder, sampleid = sampleid))
    os.system("mkdir -p {output}/HC/{sampleid}".format(output = args.outputfolder, sampleid = sampleid))
    os.system("mkdir -p {output}/logs".format(output = args.outputfolder))
    os.system("mkdir -p {output}/igv_tracks".format(output = args.outputfolder))
    os.system("mv {output}/{sampleid}/*.xml {output}/".format(sampleid = sampleid, output = args.outputfolder)) 
    os.system("mv {output}/{sampleid}/*.log {output}/logs/".format(sampleid = sampleid, output = args.outputfolder))
    os.system("mv {output}/{sampleid}/*.igv {output}/igv_tracks/".format(sampleid = sampleid, output = args.outputfolder))
    os.system("mv {output}/{sampleid}/HC*.vcf {output}/HC/".format(sampleid = sampleid, output = args.outputfolder))
    os.system("mv {output}/{sampleid}/UMCU*.vcf {output}/UMCU/".format(sampleid = sampleid, output = args.outputfolder))
    os.system("mv {output}/{sampleid}/HC* {output}/HC/{sampleid}".format(sampleid = sampleid, output = args.outputfolder))
    os.system("mv {output}/{sampleid}/UMCU* {output}/UMCU/{sampleid}".format(sampleid = sampleid, output = args.outputfolder))
    os.system("rm -r {output}/{sample}".format(output = args.outputfolder, sample = sampleid))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Path to input folder containing BAM files to process')
    parser.add_argument('outputfolder', help='Path to output folder')
    parser.add_argument('simjobs', help='number of simultanious samples to proces. Note: make sure similar threads are reseved in session!')
    parser.add_argument('--pipeline', default='iap', choices=['nf', 'iap'], help='pipeline used for sample processing (nf = nexflow, IAP = illumina analysis pipeline')
    parser.add_argument('--refset', default = settings.refset, help='Reference set to be used')
    parser.add_argument('--exomedepth', default = "/hpc/diaggen/software/production/Dx_resources/ExomeDepth/run_ExomeDepth.py", help='Full path to exomedepth script')
    parser.add_argument('--expectedCNVlength',default=settings.expectedCNVlength, help='expected CNV length (basepairs) taken into account by ExomeDepth [default expectedCNVlength in settings.py]')
    args = parser.parse_args()


    if args.pipeline == "iap":
        bams = subprocess.getoutput("find -L {0} -type f -name \"*.realigned.bam\" -not -ipath \"*exomedepth*\"".format(args.inputfolder)).split()
    elif args.pipeline == "nf":
        bams = subprocess.getoutput("find -L {0} -type f -name \"*.bam\" -not -ipath \"*exomedepth*\"".format(args.inputfolder)).split()

    print("Number of BAM files = "+str(len(bams)))
    with Pool(processes=int(args.simjobs)) as pool:
        result = pool.map(process, bams, 1)
