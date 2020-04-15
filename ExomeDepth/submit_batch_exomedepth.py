#! /usr/bin/env python
import os
from optparse import OptionParser, OptionGroup
import subprocess
import settings

if __name__ == "__main__":
    parser = OptionParser()
    group = OptionGroup(parser, "Main options")
    group.add_option("-i", dest = "input_folder", metavar = "[PATH]",
                     help = "input folder for BAM files"
                     )
    group.add_option("-o", dest = "output_folder", metavar = "[STRING]",
                     help = "Path to output folder"
                     )
    group.add_option("-e", default = "/hpc/diaggen/software/development/Dx_resources_ED_NF_python3/ExomeDepth/run_ExomeDepth.py",
                     dest = "exomedepth", metavar = "[PATH]", help = "full path to exomedepth script [default = /hpc/diaggen/software/production/Dx_resources/ExomeDepth/run_ExomeDepth.py]"
                     )
    group.add_option("-m", dest = "email", metavar = "[STRING]",
                     help = "email adress"
                     )
    #group.add_option("-s", default = " -cwd -l h_rt=4:0:0 -l h_vmem=20G -q all.q",
    #                 dest = "qsub_settings", metavar = "[STRING]", help = "qsub setting [default = qsub -cwd -l h_rt=4:0:0 -pe threaded 4 -l h_vmem=20G -q all.q]"
    #                 )
    group.add_option("-r", dest = "refset", metavar = "[STRING]",
                     help = "reference set to be used [default = reference set in setting.py]")

    parser.add_option_group(group)
    (opt, args) = parser.parse_args()

    if opt.refset:
        refset = opt.refset
    else:
        refset = settings.refset
   
    bams = subprocess.getoutput("find -L {0} -iname \"*realigned.bam\" ".format(opt.input_folder)).split()
    for bam in bams:
        #sampleid = bam.split("/")[-1].split("_")[0]
        os.system("python {exomedepth} -c -m {email} --ib={bam} -o {output} --refset={refset}".format(
            exomedepth = opt.exomedepth, 
            email = opt.email, 
            bam = bam, 
            output = opt.output_folder,
            refset = refset,
            #settings = opt.qsub_settings, 
            #sample = sampleid
        ))
