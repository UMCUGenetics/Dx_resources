#! /usr/bin/env python
import glob 
import argparse
import pysam

def detect_merge(arguments):
    write_file = open(arguments.outputfile, "w")
    bams = glob.glob("{}*/**/*.bam".format(arguments.inputfolder), recursive=True)
    for bam in bams:
        if "tmp" not in bam and "ExomeDepth" not in bam and "exomedepth" not in bam and "work" not in bam: # Skip files in folder with keywords
            workfile = pysam.AlignmentFile(bam, "rb")

            """ Extract sample_ID """
            sample_id = bam.split("/")[-1].split("_")[0]
            if "bam" in sample_id: # in case new bam file naming
                sample_id = sample_id.split(".")[0]

            """ Extract run ID """
            run_id = bam.split("/")[-3]
            if run_id == sample_id:  # Run is old IAP run
                run_id = bam.split("/")[-4]

            """ Removes deck (A/B) from run_barcode ID. This information is not present in BAM file """
            run_barcode = run_id.split("_")[3][1:]

            """ Loop reads in BAM (for maximum max_reads) and check barcode ID in BAM with run barcode ID"""
            itteration = 0
            barcode_list = {}
            for line in workfile:
                barcode = line.query_name.split(":")[2]
                if barcode not in barcode_list:
                    barcode_list[barcode] = 0
                barcode_list[barcode] += 1 
                itteration += 1
                if itteration == arguments.max_reads:
                    break
            if len(barcode_list) > 1:  # Sample consisting of multiple sequence runs is considered merge sample
                write_file.write("{sample_id}\t{run_id}\t{reason}\n".format(sample_id=sample_id,run_id=run_id,reason="multiple_runs"))
            else: 
                for item in barcode_list:
                    sample_barcode = item
                if run_barcode in barcode_list: # No merge
                    pass
                else:  # Considered merge sample
                    write_file.write("{sample_id}\t{run_id}\t{reason}\t{sample_barcode}\n".format(sample_id=sample_id,run_id=run_id,reason="different_barcode",sample_barcode=sample_barcode))
    write_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='input folder which included BAM files')
    parser.add_argument('outputfile', help='output filename of identified merge samples')
    parser.add_argument('--max_reads', default=1000, help='Maximum reads to consider in the BAM file (will speed up calculations if lower)(default = 1000)') 
    arguments = parser.parse_args()
    detect_merge(arguments)
