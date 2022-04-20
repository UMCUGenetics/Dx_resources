#! /usr/bin/env python3
import argparse
import glob
import pysam
import settings

def is_valid_read(read):
    """Check if a read is properly mapped."""
    if (read.mapping_quality >= 20 and read.reference_end and read.reference_start):
        return True
    return False

def get_gender(bam, arguments):
    """Determine chrX/chrY ratio based on read count in bam (excl PAR)."""
    with pysam.AlignmentFile(bam, "rb") as bam_file:
        y_reads = float(sum([is_valid_read(read) for read in bam_file.fetch(region=arguments.locus_y)]))
        x_reads = float(sum([is_valid_read(read) for read in bam_file.fetch(region=arguments.locus_x)]))
        total_reads = float(bam_file.mapped)
        y_ratio_perc = (y_reads / total_reads) * 100
        x_ratio_perc = (x_reads / total_reads) * 100
    
        if y_ratio_perc <= arguments.ratio_y_female and x_ratio_perc >= arguments.ratio_x_female:
            return "female\t{y_ratio_perc:.4f}\t{x_ratio_perc:.4f}".format(y_ratio_perc=y_ratio_perc, x_ratio_perc=x_ratio_perc)

        elif y_ratio_perc >= arguments.ratio_y_male and x_ratio_perc <= arguments.ratio_x_male:
            return "male\t{y_ratio_perc:.4f}\t{x_ratio_perc:.4f}".format(y_ratio_perc=y_ratio_perc, x_ratio_perc=x_ratio_perc)

        else:
            return "unknown or other sex chromosome combination\t{y_ratio_perc:.4f}\t{x_ratio_perc:.4f}".format(y_ratio_perc=y_ratio_perc, x_ratio_perc=x_ratio_perc)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Path to root folder of analysis')
    parser.add_argument('--locus_y', default=settings.locus_y, help='Coordinates for includes region on chromosome X (default = locus_y in settings.py)')
    parser.add_argument('--locus_x', default=settings.locus_x, help='Threshold for maximum allowed CNV calls (default = locus_x in settings.py)')
    parser.add_argument('--ratio_y_female', default=settings.ratio_y[0], type=float, help='Maximum Y ratio threshold females (default = ratio_y[0] in settings.py)')
    parser.add_argument('--ratio_y_male', default=settings.ratio_y[1], type=float, help='Minimum Y ratio threshold males (default = ratio_y[1] in settings.py)')
    parser.add_argument('--ratio_x_male', default=settings.ratio_x[0], type=float, help='Maximum X ratio threshold males (default = ratio_x[0] in settings.py)')
    parser.add_argument('--ratio_x_female', default=settings.ratio_x[1], type=float, help='Minimum X ratio threshold females (default = ratio_x[1] in settings.py)')
    arguments = parser.parse_args()

    bams = glob.glob("{}/**/*.bam".format(arguments.inputfolder), recursive=True)
    for bam in bams:
        print("{bam}\t{gender_result}".format(bam=bam, gender_result=get_gender(bam, arguments)))
