#! /usr/bin/env python3

import statistics
import re
import glob
import pysam


def exomedepth_summary(exomedepth_logs, stdout=None):
    stats_dic = {"CR": [], "PD": [], "TC": []}
    sample_lines = ""
    for exomedepth_qc_file in exomedepth_logs:
        with open(exomedepth_qc_file, 'r') as exomedepth_qc_lines:
            for line in exomedepth_qc_lines:
                splitline = line.split()
                correlation = float(splitline[3])
                deldupratio = float(splitline[4])
                totalcount = float(splitline[5])
                warnings = splitline[6:]
                stats_dic["CR"].append(correlation)
                stats_dic["PD"].append(deldupratio)
                stats_dic["TC"].append(totalcount)
                sample_line = (
                    "{sample};CM={model};REFSET={refset};CR={correl};PD={deldupratio};TC={totalcount}\t{warnings}\r"
                ).format(
                    sample=splitline[0],
                    model=splitline[1],
                    refset=splitline[2],
                    correl="%.4f" % correlation,
                    deldupratio="%.2f" % deldupratio,
                    totalcount="%.0f" % totalcount,
                    warnings="\t".join(warnings)
                )
                if stdout:
                    print(sample_line)
                else:
                    sample_lines = "{}\n{}".format(sample_lines, sample_line)

    mean_CR_line = "#Average_CR={}\r".format("%.4f" % statistics.mean(stats_dic["CR"]))
    mean_PD_line = "#Average_PD={}\r".format("%.2f" % statistics.mean(stats_dic["PD"]))
    mean_TC_line = "#Average_TC={}\r".format("%.2f" % statistics.mean(stats_dic["TC"]))
    median_CR = "#Median_CR={}\r".format("%.4f" % statistics.median(stats_dic["CR"]))
    median_PD = "#Median_PD={}\r".format("%.2f" % statistics.median(stats_dic["PD"]))
    median_TC = "#Median_TC={}\r".format("%.2f" % statistics.median(stats_dic["TC"]))

    if stdout:
        print((
            "\r{}\n{}\n{}\n{}\n\r\n{}\n{}\n{}"
        ).format(
            sample_lines, mean_CR_line, mean_PD_line, mean_TC_line, median_CR, median_PD, median_TC
        ))
    else:
        return (
            "\r{}\n{}\n{}\n{}\n\r\n{}\n{}\n{}"
        ).format(
            sample_lines, mean_CR_line, mean_PD_line, mean_TC_line, median_CR, median_PD, median_TC
        )


def detect_merge(inputfolder, outputfile):
    merge_file = open(outputfile, "w")
    bams = glob.glob("{}*/**/*.bam".format(inputfolder), recursive=True)
    for bam in bams:
        skip_file_keywords = ['tmp', 'exomedepth', 'work']
        if not any(keyword in bam.lower() for keyword in skip_file_keywords):  # Skip files in folder with keywords
            with pysam.AlignmentFile(bam, "rb") as bam_file:
                """ Extract sample_ID """
                sample_id = re.split('_|\.', bam.split("/")[-1])[0]
                """ Extract run ID """
                run_id = bam.split("/")[-3]
                if run_id == sample_id:  # Run is old IAP run
                    run_id = bam.split("/")[-4]

                """ Removes deck (A/B) from run_barcode ID. This information is not present in BAM file """
                run_barcode = re.sub('^[A|B]', '', run_id.split("_")[3])
                """ Compare barcode ID in BAM with run barcode ID(s) of sample """
                sample_barcodes = list(set([read_group['PU'] for read_group in bam_file.header['RG']]))
                sample_barcode = re.sub('^[A|B]', '', sample_barcodes[0])

                if len(sample_barcodes) > 1:  # Sample consisting of multiple sequence runs is considered merge sample
                    merge_file.write("{sample_id}\t{run_id}\t{reason}\n".format(
                        sample_id=sample_id, run_id=run_id, reason="multiple_runs"
                        ))
                elif run_barcode != sample_barcode:
                    merge_file.write("{sample_id}\t{run_id}\t{reason}\t{sample_barcode}\n".format(
                        sample_id=sample_id, run_id=run_id, reason="different_barcode", sample_barcode=sample_barcode
                        ))

    merge_file.close()
