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
    """Determine chrY ratio based on read count in bam (excl PAR)."""
    with pysam.AlignmentFile(bam, "rb") as bam_file:
        y_reads = float(sum([is_valid_read(read) for read in bam_file.fetch(region=arguments.gender_determination_locus_y)]))
        x_reads = float(sum([is_valid_read(read) for read in bam_file.fetch(region=arguments.gender_determination_locus_x)]))
        total_reads = float(bam_file.mapped)
        y_ratio = (y_reads / total_reads) * 100
        x_ratio = (x_reads / total_reads) * 100

        if (y_ratio <= arguments.gender_determination_y_ratio_female
           and x_ratio >= arguments.gender_determination_x_ratio_female):
            return "female\t{y_ratio:.2f}\t{x_ratio:.2f}".format(y_ratio=y_ratio, x_ratio=x_ratio)

        elif y_ratio >= arguments.gender_determination_y_ratio_male and x_ratio <= arguments.gender_determination_x_ratio_male:
            return "male\t{y_ratio:.2f}\t{x_ratio:.2f}".format(y_ratio=y_ratio, x_ratio=x_ratio)

        else:
            return "unknown or other sex chromosome combination\t{y_ratio}\t{x_ratio}".format(y_ratio=y_ratio, x_ratio=x_ratio)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Path to root folder of analysis')
    parser.add_argument(
        '--gender_determination_locus_y',
        default=settings.gender_determination_locus_y,
        help='Coordinates for includes region on chromosome X (default = gender_determination_locus_y in settings.py)'
    )
    parser.add_argument(
        '--gender_determination_locus_x',
        default=settings.gender_determination_locus_x,
        help='Threshold for maximum allowed CNV calls (default = gender_determination_locus_x in settings.py)'
    )
    parser.add_argument(
        '--gender_determination_y_ratio_female',
        default=settings.gender_determination_y_ratio[0],
        type=float,
        help='Maximum Y ratio threshold females (default = gender_determination_y_ratio[0] in settings.py)'
    )
    parser.add_argument(
        '--gender_determination_y_ratio_male',
        default=settings.gender_determination_y_ratio[1],
        type=float,
        help='Minimum Y ratio threshold males (default = gender_determination_y_ratio[1] in settings.py)'
    )
    parser.add_argument(
        '--gender_determination_x_ratio_male',
        default=settings.gender_determination_x_ratio[0],
        type=float,
        help='Maximum X ratio threshold males (default = gender_determination_x_ratio[0] in settings.py)'
    )
    parser.add_argument(
        '--gender_determination_x_ratio_female',
        default=settings.gender_determination_x_ratio[1],
        type=float,
        help='Minimum X ratio threshold females (default = gender_determination_x_ratio[1] in settings.py)'
    )
    arguments = parser.parse_args()

    bams = glob.glob("{}/**/*.bam".format(arguments.inputfolder), recursive=True)
    for bam in bams:
        print("{bam}\t{gender_result}".format(bam=bam, gender_result=get_gender(bam, arguments)))
