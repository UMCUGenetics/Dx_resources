#! /usr/bin/env python
import re
import glob 
import argparse
import pysam

def detect_merge(arguments):
    write_file = open(arguments.outputfile, "w")
    bams = glob.glob("{}*/**/*.bam".format(arguments.inputfolder), recursive=True)
    for bam in bams:
        skip_file_keywords = ['tmp', 'exomedepth', 'work']
        if not any(keyword in bam.lower() for keyword in skip_file_keywords): # Skip files in folder with keywords
            bam_file = pysam.AlignmentFile(bam, "rb")

            """ Extract sample_ID """
            sample_id = re.split('_|\.', bam.split("/")[-1])[0]

            """ Extract run ID """
            run_id = bam.split("/")[-3]
            if run_id == sample_id:  # Run is old IAP run
                run_id = bam.split("/")[-4] 

            """ Removes deck (A/B) from run_barcode ID. This information is not present in BAM file """
            run_barcode = run_id.split("_")[3][1:]

            """ Compare barcode ID inin BAM with run barcode ID(s) of sample """
            sample_barcodes = list(set([read_group['PU'] for read_group in bam_file.header['RG']]))
            sample_barcode = re.sub('^[A|B]', '', sample_barcodes[0])

            if len(sample_barcodes) > 1:  # Sample consisting of multiple sequence runs is considered merge sample
                write_file.write("{sample_id}\t{run_id}\t{reason}\n".format(
                sample_id=sample_id, run_id=run_id, reason="multiple_runs"
                ))
            elif run_barcode != sample_barcode:
                write_file.write("{sample_id}\t{run_id}\t{reason}\t{sample_barcode}\n".format(
                sample_id=sample_id, run_id=run_id, reason="different_barcode", sample_barcode=sample_barcode
                ))

    write_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='input folder which included BAM files')
    parser.add_argument('outputfile', help='output filename of identified merge samples')
    arguments = parser.parse_args()
    detect_merge(arguments)
