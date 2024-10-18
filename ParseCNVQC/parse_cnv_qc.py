#! /usr/bin/env python3

import argparse
import os
import time
from email import encoders
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
import smtplib
import mimetypes
from datetime import date

import settings


def send_email(sender, receivers, subject, text, attachment=None):
    """Send email utility function."""
    mail = MIMEMultipart()
    mail['Subject'] = subject
    mail['From'] = sender
    mail['To'] = ';'.join(receivers)
    if attachment:
        filename = attachment.split('/')[-1]
        fp = open(attachment, 'rb')
        ctype, encoding = mimetypes.guess_type(attachment)
        if ctype is None or encoding is not None:
            """ No guess could be made, or the file is encoded (compressed), so use a generic bag-of-bits type."""
            ctype = 'application/octet-stream'
        maintype, subtype = ctype.split('/', 1)
        msg = MIMEBase(maintype, subtype)
        msg.set_payload(fp.read())
        fp.close()
        """ Encode the payload using Base64 """
        encoders.encode_base64(msg)
        msg.add_header('Content-Disposition', 'attachment', filename=filename)
        mail.attach(msg)
    msg = MIMEText(text, 'html')
    mail.attach(msg)
    m = smtplib.SMTP('smtp-open.umcutrecht.nl')
    m.sendmail(sender, receivers, mail.as_string())
    m.quit()


def make_mail(today, daysago, attachment, run_status):
    subject = 'CNV QC summary {}'.format(today)
    text = '<html><body><p> CNV QC summary {} from the last {} days'.format(today, daysago)
    text = '{}<br><br>Likely complete runs:'.format(text)
    if run_status["full"]:
        for run in run_status["full"]:
            text = '{}<br>&nbsp;&nbsp;&nbsp;&nbsp;{}'.format(text, run)
    else:
        text = '{}<br>&nbsp;&nbsp;&nbsp;&nbsp;none'.format(text)
    text = '{}<br><br>Likely partial runs:'.format(text)
    if run_status["partial"]:
        for run in run_status["partial"]:
            text = '{}<br>&nbsp;&nbsp;&nbsp;&nbsp;{}'.format(text, run)
    else:
        text = '{}<br>&nbsp;&nbsp;&nbsp;&nbsp;none'.format(text)
    text = '{}<br><br><b>This analysis is not for diagnostic use.</b>'.format(text)
    text = '{}</p></body></html>'.format(text)
    send_email(settings.email_from, settings.email_to, subject, text, attachment)


def include_sample_number(folder, rawfolder, projects, warnings):
    """
    Determine number of samples in a run based on the SampleSheet.
    Only includes samples that are in predefined projects as stated in settings file.
    Includes a warning if Samplesheet.csv is not detected in the raw data folder

    Args:
        folder (string): full path to raw data folder
        rawfolder (dict):
            key: full path to raw data folder, values: list of all processed folder from run, total number of samples (int)
        projects (list): projectIDs to include in total sample calculation. Excludes i.e. non-dx projects in calculations
        warnings (list): List with warning messages

    Returns:
        rawfolder (dict): rawfolder including sample count
        warnings (list): message list including potential new warnings
    """
    number_samples_run = 0
    lanes = []
    lane_index = ""
    if os.path.exists("{}/SampleSheet.csv".format(folder)):
        with open("{}/SampleSheet.csv".format(folder), 'r') as samplesheet:
            sample_section = False
            for line in samplesheet:
                if sample_section:
                    for project in projects:
                        if project in line.upper():
                            number_samples_run += 1
                            if line.split(",")[lane_index] not in lanes:
                                lanes.append(line.split(",")[lane_index])
                if "Sample_ID" not in line:
                    continue
                else:
                    sample_section = True
                    header = [column for column in line.split(",")]
                    lane_index = header.index('Lane')
    else:
        warnings.append("no samplesheet for run {}, assuming unknown number of samples in run".format(folder))

    # prevent division by zero.
    if len(lanes) > 0:
        rawfolder[folder][1] += number_samples_run/len(lanes)
    return rawfolder, warnings


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--daysago',
        default=settings.daysold,
        type=int,
        help='number of days ago to be included (default = {})'.format(settings.daysold)
    )
    parser.add_argument(
        '--outputfolder',
        default=settings.outputfolder,
        help='output folder of file (default = {})'.format(settings.outputfolder)
    )
    parser.add_argument(
        '--inputfolder',
        default=settings.inputfolder,
        help='path to processed data folder (default = {})'.format(settings.inputfolder)
    )
    parser.add_argument(
        '--rawsdatafolder',
        default=settings.rawdatafolder,
        help='path to raw data folder (default = {})'.format(settings.rawdatafolder)
    )

    parser.add_argument(
        '--threshold',
        default=settings.threshold,
        type=float,
        help='threshold to include warning (#assesed < #samplesheet (default = {})'.format(settings.threshold)
    )

    parser.add_argument(
        '--windows_nl',
        action='store_true',
        help='make and send Windows/Dutch specific file (do changed to comma). [Default = off]'
    )

    args = parser.parse_args()

    runs = {}
    warnings = []
    folders = sorted([file for file in os.listdir(args.inputfolder) if os.path.isdir(os.path.join(args.inputfolder, file))])
    days_ago = time.time() - args.daysago * 24 * 60 * 60
    for folder in folders:
        folder_timestamp = os.stat("{}/{}".format(args.inputfolder, folder)).st_mtime
        if folder_timestamp > days_ago:
            run = "_".join(folder.split("_")[0:4])
            if run not in runs:
                runs[run] = ''

    rawfolder = {}
    for folder in folders:
        run = "_".join(folder.split("_")[0:4])
        if run in runs and "old" not in folder:
            run_path = "{}/{}".format(args.rawsdatafolder, run)
            if run_path not in rawfolder:
                rawfolder[run_path] = [[], 0]
            rawfolder[run_path][0].append(folder)

    # Get number of samples for each run
    for folder in rawfolder:
        rawfolder, warnings = include_sample_number(folder, rawfolder, settings.projects, warnings)

    folder_summary = {}
    sample_qc = []
    for folder in rawfolder:
        total_samples = rawfolder[folder][1]
        count = {"total": 0, "CR": 0, "TC": 0, "PD": 0, "QC_failed": 0}
        for analysis in rawfolder[folder][0]:
            qc_file = "{path}/{analysis}/QC/CNV/{analysis}_exomedepth_summary.txt".format(
                path=args.inputfolder, analysis=analysis
            )
            if os.path.exists(qc_file):
                with open(qc_file, 'r') as qc_lines:
                    for line in qc_lines:
                        if "REFSET" in line:
                            splitline = line.rstrip().split("\t")[0].split(";")
                            count["total"] += 1
                            CR_passed = "Failed"
                            PD_passed = "Failed"
                            TC_passed = "Failed"
                            QC_failed = "Ok"
                            for item in splitline:
                                splititem = item.split("=")
                                if splititem[0] in settings.QC:
                                    if splititem[0] == "CR" and float(splititem[1]) < settings.QC["CR"][0]:
                                        count["CR"] += 1
                                        CR = float(splititem[1])
                                        QC_failed = "Failed"
                                    elif splititem[0] == "CR":
                                        CR_passed = "Ok"
                                        CR = float(splititem[1])
                                    if (splititem[0] == "PD"
                                       and (float(splititem[1]) < settings.QC["PD"][0]
                                       or float(splititem[1]) > settings.QC["PD"][1])):
                                        count["PD"] += 1
                                        PD = float(splititem[1])
                                        QC_failed = "Failed"
                                    elif splititem[0] == "PD":
                                        PD_passed = "Ok"
                                        PD = float(splititem[1])
                                    if (splititem[0] == "TC"
                                       and (float(splititem[1]) < settings.QC["TC"][0]
                                       or float(splititem[1]) > settings.QC["TC"][1])):
                                        count["TC"] += 1
                                        TC = float(splititem[1])
                                        QC_failed = "Failed"
                                    elif splititem[0] == "TC":
                                        TC_passed = "Ok"
                                        TC = float(splititem[1])
                            if QC_failed == "Failed":
                                count["QC_failed"] += 1
                            sample_qc.append("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                                analysis, "\t".join(splitline[0:3]), CR, PD, TC, CR_passed, PD_passed, TC_passed, QC_failed)
                            )

            else:
                warnings.append("no QC CNV summary file is detected for analysis {}, skipping".format(analysis))

        if folder not in folder_summary:
            folder_summary[folder] = [count, total_samples]

    today = date.today().strftime("%y_%m_%d")
    output = "{}/CNV_QC_summary_{}.tsv".format(args.outputfolder, today)
    run_status = {"full": [], "partial": []}
    with open(output, 'w') as outputfile:
        outputfile.write("Note: The software that produced this file has not been validated for diagnostic use.\n\n")
        outputfile.write("###SAMPLE specific QC\n")
        outputfile.write("Analysis\tSample\tHC\tREFSET\tCR\tPD\tTC\tCR_passed\tPD_passed\tTC_passed\tQC_passed\n")
        for sample in sample_qc:
            outputfile.write("{}\n".format(sample))
        outputfile.write("\n###Run specific QC\n")
        outputfile.write("Run\tFailed%\tPassed%\t#SamplesProject\t#SamplesRun\n")
        for folder in folder_summary:
            total_assesed = float(folder_summary[folder][0]["total"])
            total_run = float(folder_summary[folder][1])
            total_failed = float(folder_summary[folder][0]["QC_failed"])

            failed_perc = float("%.2f" % ((total_failed/total_assesed) * 100))
            passed_perc = 100 - failed_perc
            outputfile.write("{}\t{}\t{}\t{}\t{}".format(
                folder.split("/")[-1], failed_perc, passed_perc, int(total_assesed), int(total_run)
            ))
            if total_run > 0:
                if (total_assesed/total_run) < args.threshold:
                    outputfile.write((
                        "\tWARNING: number of samples within run is less than {}% of total in samplesheet."
                    ).format(
                        args.threshold * 100
                    ))
                    run_status["partial"].append(folder.split("/")[-1])
                else:
                    run_status["full"].append(folder.split("/")[-1])
            else:
                outputfile.write("\tWARNING: unknown number of samples in run")
                run_status["partial"].append(folder.split("/")[-1])
            outputfile.write("\n")

        outputfile.write("\n###WARNINGS\n")
        for warning in warnings:
            outputfile.write("WARNING:{}\n".format(warning))

    if args.windows_nl:
        output_windows = "{}/CNV_QC_summary_{}_win.tsv".format(args.outputfolder, today)
        os.system("sed 's/\./,/g' {} > {}".format(output, output_windows))
        make_mail(today, args.daysago, output_windows, run_status)
    else:
        make_mail(today, args.daysago, output, run_status)
