#! /usr/bin/env python
import sys, os, re
import commands
import getopt

## Log git commits of used repo's
os.system(" date >> GIT_IAP.log")
os.system(" date >> GIT_ExonCov.log")
os.system(" date >> GIT_Dx_resources.log")
os.system("git --git-dir=/hpc/cog_bioinf/common_scripts/IAP_Dx_TARGET/.git log > GIT_IAP.log")
os.system("git --git-dir=/hpc/cog_bioinf/common_scripts/ExonCov/.git log > GIT_ExonCov.log")
os.system("git --git-dir=/hpc/cog_bioinf/data/mapping/diagnostiek/Dx_resources/.git log > GIT_Dx_resources.log")

# Make run_stats.txt file from all flagstats
os.system("/hpc/cog_bioinf/data/mapping/diagnostiek/Dx_resources/get_stats_from_flagstat.pl >run_stats.txt")

sys.exit("Finished")

