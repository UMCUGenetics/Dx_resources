#! /usr/bin/env python3
import sys
import os
from statistics import mean
from statistics import stdev

wkdir = sys.argv[1]  #Input folder

""" Make bam_coverage count for each BAM file first"""
os.system("paste {0}/*bam_coverage > {0}/full_table.txt".format(wkdir))  #Make one table of all bam_coverage counts

with open(str(wkdir) + "/full_table.txt", "r") as infile:
    for line in infile: 
        if "#" not in line:
            splitline = line.split("\t") 
            read_list = []
            depth_list = []

            """Read depth stats."""
            mean_coverage_column = 4
            while mean_coverage_column < len(splitline):
                depth_list += [float(splitline[mean_coverage_column])]
                mean_coverage_column += 6  #Go to next meancoverage column. note: this dependants on the sambamaba output file!

            """Calculate mean, std, cv."""
            mean_depth = float(mean(depth_list))
            std_depth = float(stdev(depth_list))
            if mean_depth > 0:
                cv_depth = (std_depth/mean_depth) * 100
            else:
                cv_depth = "inf"
      
            """Print only targets that are within requirements"""
            print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(
               splitline[0].rstrip(),
               splitline[1].rstrip(),
               splitline[2].rstrip(),
               "%.2f" % float(mean_depth),
               "%.2f" % float(std_depth),
               "%.2f" % float(cv_depth)
            ))

