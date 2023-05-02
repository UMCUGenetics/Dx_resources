#! /usr/bin/env python3

import os
import argparse
import subprocess
import glob
import sys
from multiprocessing import Pool
from datetime import date

import database.functions
import utils.utils
import settings


def exomedepth_analysis(bam, args, gender_dic, suffix_dic, refset_dic):
    bamfile = bam.split("/")[-1]
    sampleid = database.functions.get_sample_id(bam)
    bam_path = os.path.abspath(bam)
    run = args.inputfolder.rstrip("/").split("/")[-1]

    if args.reanalysis and sampleid in refset_dic:  # Use refset in reanalysis argument file if provided in argument
        refset = refset_dic[sampleid]
    else:
        refset = database.functions.add_sample_to_db_and_return_refset_bam(bam_path, settings.refset)[2]

    os.makedirs(f"{args.outputfolder}/{sampleid}", exist_ok=True)
    os.chdir(f"{args.outputfolder}/{sampleid}")
    os.symlink(f"{bam}", f"{args.outputfolder}/{sampleid}/{bamfile}")
    os.symlink(f"{bam}.bai", f"{args.outputfolder}/{sampleid}/{bamfile}.bai")

    action = (
        f"python {args.exomedepth} callcnv {args.outputfolder}/{sampleid} {bamfile} {run} {sampleid} "
        "--refset {refset} --expectedCNVlength {args.expectedCNVlength} --pipeline {args.pipeline}"
        )

    if sampleid in gender_dic:
        gender = gender_dic[sampleid]
        action = f"{action} --refset_gender {gender}"

    if sampleid in suffix_dic:
        suffix = suffix_dic[sampleid]
        action = f"{action} --vcf_filename_suffix {suffix}"

    os.system(action)

    os.chdir(f"{args.outputfolder}")

    os.makedirs(f"{args.outputfolderoutput}/logs", exist_ok=True)
    os.makedirs(f"{args.outputfolder}/igv_tracks", exist_ok=True)
    os.makedirs(f"{args.outputfolder}/UMCU/", exist_ok=True)
    os.makedirs(f"{args.outputfolder}/HC/", exist_ok=True)

    os.system(f"mv {args.outputfolder}/{sampleid}/*.log {args.outputfolder}/logs/")
    os.system(f"mv {args.outputfolder}/{sampleid}/*.igv {args.outputfolder}/igv_tracks/")
    os.system(f"mv {args.outputfolder}/{sampleid}/HC*.vcf {args.outputfolder}/HC/")
    os.system(f"mv {args.outputfolder}/{sampleid}/UMCU*.vcf {args.outputfolder}/UMCU/")
    os.system(f"rm -r {args.outputfolder}/{sampleid}")

    hc_cnv_vcf = f"{args.outputfolder}/HC/HC_{refset}_{sampleid}_{run}_exome_calls"
    ed_igv = f"{args.outputfolder}/igv_tracks/HC_{refset}_{sampleid}_{run}"
    if sampleid in suffix_dic:
        hc_cnv_vcf = f"{hc_cnv_vcf}_{suffix_dic[sampleid]}.vcf"
        ed_igv = f"{ed_igv}_{suffix_dic[sampleid]}_ref.igv"
    else:
        hc_cnv_vcf = f"{hc_cnv_vcf}.vcf"
        ed_igv = f"{ed_igv}_ref.igv"

    return [sampleid, bamfile, refset, hc_cnv_vcf, ed_igv]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Path to root folder of analysis')
    parser.add_argument('outputfolder', help='Path to output folder')
    parser.add_argument('runid', help='Run ID')
    parser.add_argument(
        'simjobs',
        help='number of simultaneous samples to proces. Note: make sure similar threads are reserved in session!'
    )
    parser.add_argument('pedfile', help='full path to ped file')
    parser.add_argument(
        '--pipeline', default='nf', choices=['nf', 'iap'],
        help='pipeline used for sample processing (nf = nexflow (default), IAP = illumina analysis pipeline)'
    )
    parser.add_argument(
        '--reanalysis',
        help='Tab delimited file with SampleID, RefsetID, and optional reanalysis female/male mode (see settings.py)'
        'to be used. If samples are not present in the file, default refset and female/male is used'
    )
    parser.add_argument(
        '--exomedepth',
        default="{}/run_ExomeDepth.py".format(settings.cwd), help='Full path to exomedepth script'
    )
    parser.add_argument(
        '--expectedCNVlength', default=settings.expectedCNVlength,
        help='expected CNV length (basepairs) taken into account by ExomeDepth [default expectedCNVlength in settings.py]'
    )
    args = parser.parse_args()

    today = date.today().strftime("%d%m%y")
    user = subprocess.getoutput("whoami")

    """Find BAM files to be processed"""
    if args.pipeline == "iap":
        bam_files = (
            set(glob.glob("{}/**/*.realigned.bam".format(args.inputfolder), recursive=True))
            - set(glob.glob("{}/[eE]xome[dD]epth*/**/*.realigned.bam".format(args.inputfolder), recursive=True))
        )
    elif args.pipeline == "nf":
        bam_files = glob.glob("{}/bam_files/**/*.bam".format(args.inputfolder), recursive=True)

    print("Number of BAM files = "+str(len(bam_files)))

    if bam_files:
        os.system("echo \"{user} {today}\tExomeDepth reanalysis performed\" >> {inputfolder}/logbook.txt".format(
            user=user, today=today, inputfolder=args.inputfolder)
        )
    else:
        sys.exit("no bam files detected")

    snv_vcf_files = []
    for snv_vcf_file in glob.glob(f"{args.inputfolder}/single_sample_vcf/*.vcf", recursive=True):
        snv_vcf_files.append(os.path.basename(snv_vcf_file))

    baf_files = []
    for baf_file in glob.glob(f"{args.inputfolder}/baf/*.igv", recursive=True):
        baf_files.append(os.path.basename(baf_file))

    upd_files = []
    for upd_file in glob.glob(f"{args.inputfolder}/upd/*.igv", recursive=True):
        upd_files.append(os.path.basename(upd_file))

    if any(len(file_list) != len(bam_files) for file_list in [snv_vcf_files, baf_files]):
        sys.exit((
            "unequal number of files found: bams={0}, snvVCF={1}, baf={2}"
        ).format(
            len(bam_files), len(snv_vcf_files), len(baf_files)
        ))

    """ Check if previous exomedepth analysis is already present """
    if os.path.isdir(args.outputfolder):
        print("Run folder not empty: assuming exomedepth has been runned. Current exomedepth folder will be archived")
        """Check if archive folder exists. If this is the case, relative path in IGV session should not be changed."""
        archivefolder = False
        archive_folder_exomedepth = "{outputfolder}/archive_{today}/".format(outputfolder=args.outputfolder, today=today)
        if os.path.isdir(archive_folder_exomedepth):
            archivefolder = True
        else:
            os.makedirs(archive_folder_exomedepth, exist_ok=True)

        """ Move original data to archive folder """
        os.system("mv {outputfolder}/* {archive_folder_exomedepth}".format(
            outputfolder=args.outputfolder, archive_folder_exomedepth=archive_folder_exomedepth
        ))

        """ Rename relative paths in IGV sessions"""
        if not archivefolder:
            os.system("sed -i 's/\"\.\.\//\"\.\.\/\.\.\//g' {archive_folder_exomedepth}/*xml".format(
                archive_folder_exomedepth=archive_folder_exomedepth
            ))

        """ Check if archive folder is already present in archived folder, and move the to the correct location."""
        if glob.glob("{outputfolder}/archive_{today}/archive*".format(outputfolder=args.outputfolder, today=today)):
            os.system("mv {outputfolder}/archive_{today}/archive* {outputfolder}/".format(
                outputfolder=args.outputfolder, today=today
            ))

        """ Copy CNV summary file into archive folder """
        if(glob.glob("{inputfolder}/QC/CNV/{runid}_exomedepth_summary.txt".format(
           inputfolder=args.inputfolder, runid=args.runid))):
            os.makedirs(
                "{inputfolder}/QC/CNV/archive_{today}/".format(
                    inputfolder=args.inputfolder, today=today
                ),
                exist_ok=True
            )

            os.system("mv {inputfolder}/QC/CNV/{runid}_exomedepth_summary.txt {inputfolder}/QC/CNV/archive_{today}/".format(
                inputfolder=args.inputfolder, runid=args.runid, today=today
            ))
    else:
        os.makedirs(args.outputfolder, exist_ok=True)

    refset_dic = {}
    gender_dic = {}
    suffix_dic = {}
    if args.reanalysis:  # Use predetermined refset for each sample based on the reanalysis option input
        with open(args.reanalysis) as refset_file:
            for line in refset_file:
                splitline = line.split()
                sampleid = splitline[0]
                refset = splitline[1]
                if sampleid not in refset_dic:
                    refset_dic[sampleid] = refset
                if len(splitline) > 2:
                    tag = splitline[2]
                    if settings.reanalysis_dic[tag]:
                        if sampleid not in gender_dic:
                            gender_dic[sampleid] = settings.reanalysis_dic[tag][0]
                        if sampleid not in suffix_dic:
                            suffix_dic[sampleid] = settings.reanalysis_dic[tag][1]
                    else:
                        sys.exit(
                            "Warning: reanalysis tag {0} in file {1} is unknown within settings.reanalysis_dic."
                            " Please change or add reanalysis tag".format(tag, args.reanalysis)
                        )

    """Start exomedepth re-analysis"""
    with Pool(processes=int(args.simjobs)) as pool:
        sampleinfo = pool.starmap(exomedepth_analysis, [[bam, args, gender_dic, suffix_dic, refset_dic] for bam in bam_files])

    """ Make CNV summary file """
    logs = glob.glob("{outputfolder}/logs/HC*stats.log".format(outputfolder=args.outputfolder), recursive=True)
    os.makedirs("{inputfolder}/QC/CNV/".format(inputfolder=args.inputfolder), exist_ok=True)

    write_summary = "{inputfolder}/QC/CNV/{runid}_exomedepth_summary.txt".format(
        inputfolder=args.inputfolder,
        runid=args.runid
    )

    with open(write_summary, 'w') as write_file_summary:
        write_file_summary.write(utils.utils.exomedepth_summary(logs))

    """ Make single sample IGV sessions for all samples """
    for analysis in sampleinfo:
        action = "python {0} single_igv {1} {2} {3} {4} --bam {5}".format(
            settings.igv_xml, args.outputfolder, analysis[0], args.runid, analysis[2], analysis[1]
        )
        if analysis[0] in suffix_dic:
            action = f"{action} --reanalysis {args.reanalysis}"
        os.system(action)

    """ Make family IGV session(s)"""
    cnv_vcf_files = []
    igv_files = []
    for analysis in sampleinfo:
        cnv_vcf_files.append(os.path.basename(analysis[3]))
        igv_files.append(os.path.basename(analysis[4]))

    bam_file_basenames = []
    for bam in bam_files:
        bam_file_basenames.append(os.path.basename(bam))

    action = (
        "python {0} family_igv {1} {2} {3} --bam_files {4} "
        "--snv_vcf_files {5} --cnv_vcf_files {6} --igv_files {7} "
        "--upd_files {8} --baf_files {9}"
    ).format(
        settings.igv_xml,
        args.outputfolder,
        args.pedfile,
        args.runid,
        " ".join(bam_file_basenames),
        " ".join(snv_vcf_files),
        " ".join(cnv_vcf_files),
        " ".join(igv_files),
        " ".join(upd_files),
        " ".join(baf_files)
    )

    if args.reanalysis:
        f"{action} --reanalysis {args.reanalysis}"

    os.system(action)
