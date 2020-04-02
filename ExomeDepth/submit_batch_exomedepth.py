#! /usr/bin/env python
import os
import sys
from optparse import OptionParser, OptionGroup
import commands

if __name__ == "__main__":
    parser = OptionParser()
    group = OptionGroup(parser, "Main options")
    group.add_option("-i", dest = "input_folder", metavar = "[PATH]",
                     help = "input folder for BAM files"
                     )
    group.add_option("-o", dest = "output_folder", metavar = "[STRING]",
                     help = "Path to output folder"
                     )
    group.add_option("-e", default = "/hpc/diaggen/software/production/Dx_resources/ExomeDepth/run_ExomeDepth.py",
                     dest = "exomedepth", metavar = "[PATH]", help = "full path to exomedepth script [default = /hpc/diaggen/software/production/Dx_resources/ExomeDepth/run_ExomeDepth.py]"
                     )
    group.add_option("-m", dest = "email", metavar = "[STRING]",
                     help = "email adress"
                     )
    group.add_option("-s", default = "qsub -cwd -l h_rt=4:0:0 -l h_vmem=20G -M m.elferink@umcutrecht.nl -q all.q",
                     dest = "settings", metavar = "[STRING]", help = "qsub setting [default = qsub -cwd -l h_rt=4:0:0 -pe threaded 4 -l h_vmem=20G -q all.q]"
                     )
    parser.add_option_group(group)
    (opt, args) = parser.parse_args()

    input_folder = opt.input_folder
    output_folder = opt.output_folder
    settings = opt.settings
    exomedepth = opt.exomedepth
    email = opt.email

    bams = commands.getoutput("find -L {0} -iname \"*realigned.bam\" ".format(input_folder)).split()
    for bam in bams:
        sampleid = bam.split("/")[-1].split("_")[0]
        os.system("echo \"python {0} -c -m {1} --ib={2} -o {3}\" | {4} -M {1} -N {5}.sh".format(exomedepth, email, bam, output_folder, settings, sampleid))

