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

def get_gender(bam, arguments):
    """Determine chrY ratio based on read count in bam (excl PAR)."""
    workfile = pysam.AlignmentFile(bam, "rb")
    y_reads = float(sum([valid_read(read) for read in workfile.fetch(region=arguments.locusy)]))
    x_reads = float(sum([valid_read(read) for read in workfile.fetch(region=arguments.locusx)]))
    total_reads = float(workfile.mapped)
    y_ratio = float("%.2f" % ((y_reads / total_reads) * 100))
    x_ratio = float("%.2f" % ((x_reads / total_reads) * 100))
    
    if y_ratio <= arguments.y_ratio_female and x_ratio >= arguments.x_ratio_female:
        return "female\t{y_ratio}\t{x_ratio}".format(y_ratio=y_ratio, x_ratio=x_ratio)
    elif y_ratio >= arguments.y_ratio_male and x_ratio <= arguments.x_ratio_male:
        return "male\t{y_ratio}\t{x_ratio}".format(y_ratio=y_ratio, x_ratio=x_ratio)
    else:
        return "unknown or other sex chromosome combination\t{y_ratio}\t{x_ratio}".format(y_ratio=y_ratio, x_ratio=x_ratio)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Path to root folder of analysis')
    parser.add_argument('--locusy', default=settings.locusy, help='Coordinates for includes region on chromosome X (default = locusy in settings.py)')
    parser.add_argument('--locusx', default=settings.locusx, help='Threshold for maximum allowed CNV calls (default = locusx in settings.py)')
    parser.add_argument('--y_ratio_female', default=settings.y_ratio[0], type=float, help='Maximum Y ratio threshold females (default = y_ratio[0] in settings.py)')
    parser.add_argument('--y_ratio_male', default=settings.y_ratio[1], type=float, help='Minimum Y ratio threshold males (default = y_ratio[1] in settings.py)')
    parser.add_argument('--x_ratio_male', default=settings.x_ratio[0], type=float, help='Maximum X ratio threshold males (default = x_ratio[0] in settings.py)')
    parser.add_argument('--x_ratio_female', default=settings.x_ratio[1], type=float, help='Minimum X ratio threshold females (default = x_ratio[1] in settings.py)')
    arguments = parser.parse_args()

    bams = glob.glob("{}/**/*.bam".format(arguments.inputfolder), recursive=True)
    for bam in bams:
        print("{bam}\t{gender_result}".format(bam=bam, gender_result=get_gender(bam, arguments)))
