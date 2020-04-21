#! /usr/bin/env python3

import os
import argparse
import subprocess
import settings

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Path to input folder containing BAM files to process')
    parser.add_argument('outputfolder', help='Path to output folder')
    parser.add_argument('--refset', default = settings.refset, help='Reference set to be used')
    parser.add_argument('--exomedepth', default = "/hpc/diaggen/software/production/Dx_resources/ExomeDepth/run_ExomeDepth.py", help='Full path to exomedepth script')
    args = parser.parse_args()

    bams = subprocess.getoutput("find -L {0} -iname \"*.realigned.bam\"".format(args.inputfolder)).split()
    for bam in bams:
        action = ("python {exomedepth} callcnv {output} {inputbam} {run} {sample} {refset} ".format(
            exomedepth = args.exomedepth, 
            output = args.outputfolder, 
            inputbam = bam,
            run = args.inputfolder.split("/")[-1],
            sample = bam.split("/")[-1].split("_")[0],
            refset = args.refset
        ))
        os.system(action)
