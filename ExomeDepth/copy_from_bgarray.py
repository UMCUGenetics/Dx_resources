#! /usr/bin/env python
import sys
import os
import os.path

if not os.path.isdir("./BAMS"):
    os.system("mkdir ./BAMS")

lines = open(sys.argv[1],"r").readlines() # Input file with runs to be transferred
for line in lines:
    action = "rsync -rahu --progress  --prune-empty-dirs --include=\'**/\' --include=\'**/exomedepth/**\' --include=\'**/ExomeDepth/**\' --include=\"*.bam\" --include=\"*.bai\" --include=\"logbook.txt\"  --exclude=\'*\' /data/DIT-bgarray/Illumina/Exomes/{run}  BAMS/".format(run=line.split()[0])
    printi(action)
    os.system(action)
