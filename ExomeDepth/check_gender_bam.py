#! /usr/bin/env python3
import argparse
import glob
import pysam
import settings

def valid_read(read):
    """Check if a read is properly mapped."""
    if (read.mapping_quality >= 20 and read.reference_end and read.reference_start):
        return True
    else:
        return False

def get_gender(bam):
    """Determine chrY ratio based on read count in bam (excl PAR)."""
    workfile = pysam.AlignmentFile(bam, "rb")
    locusY = settings.locusY
    locusX = settings.locusX
    yreads = float(sum([valid_read(read) for read in workfile.fetch(region=locusY)]))
    xreads = float(sum([valid_read(read) for read in workfile.fetch(region=locusX)]))
    total = float(workfile.mapped)
    yratio = float("%.2f" % ((yreads / total) * 100))
    xratio = float("%.2f" % ((xreads / total) * 100))

    if yratio <= float(settings.yratio[0]) and xratio >= float(settings.xratio[1]):
        return "female", yratio, xratio
    elif yratio >= float(settings.yratio[1]) and xratio <= float(settings.xratio[0]):
        return "male", yratio, xratio
    else:
        return "unknown or other sex chromosome combination", yratio, xratio

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Path to root folder of analysis')
    args = parser.parse_args()

    bams = glob.glob("{}/**/*.bam".format(args.inputfolder), recursive=True)
    for bam in bams:
        print(bam, *get_gender(bam))
