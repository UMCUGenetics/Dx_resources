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

    if args.reanalysis and sampleid in refset_dic:  # Use refset in reanalysis argument file if provided in argument
        refset = refset_dic[sampleid]
    else:
        refset = database.functions.add_sample_to_db_and_return_refset_bam(bam_path, settings.refset)[2]

    os.makedirs("{output}/{sample}".format(output=args.outputfolder, sample=sampleid), exist_ok=True)
    os.chdir("{output}/{sample}".format(output=args.outputfolder, sample=sampleid))

    os.symlink(
        "{bam}".format(bam=bam),
        "{output}/{sample}/{bamfile}".format(output=args.outputfolder, bamfile=bamfile, sample=sampleid)
    )

    os.symlink(
        "{bam}.bai".format(bam=bam),
        "{output}/{sample}/{bamfile}.bai".format(output=args.outputfolder, bamfile=bamfile, sample=sampleid)
    )

    action = (
        "python {exomedepth} callcnv {output}/{sample} {inputbam} {run} {sample} "
        "--refset {refset} --expectedCNVlength {length} --pipeline {pipeline}").format(
            exomedepth=args.exomedepth,
            output=args.outputfolder,
            inputbam=bamfile,
            run=args.inputfolder.rstrip("/").split("/")[-1],
            sample=sampleid,
            refset=refset,
            length=args.expectedCNVlength,
            pipeline=args.pipeline
        )

    if sampleid in gender_dic:
        gender = gender_dic[sampleid]
        action = "{action} --refset_gender {gender}".format(action=action, gender=gender)

    if sampleid in suffix_dic:
        suffix = suffix_dic[sampleid]
        action = "{action} --vcf_filename_suffix {suffix}".format(action=action, suffix=suffix)

    os.system(action)

    os.chdir("{output}".format(output=args.outputfolder))

    os.makedirs("{output}/logs".format(output=args.outputfolder), exist_ok=True)
    os.makedirs("{output}/igv_tracks".format(output=args.outputfolder), exist_ok=True)
    os.makedirs("{output}/UMCU/".format(output=args.outputfolder), exist_ok=True)
    os.makedirs("{output}/HC/".format(output=args.outputfolder), exist_ok=True)

    os.system("mv {output}/{sampleid}/*.log {output}/logs/".format(sampleid=sampleid, output=args.outputfolder))
    os.system("mv {output}/{sampleid}/*.igv {output}/igv_tracks/".format(sampleid=sampleid, output=args.outputfolder))
    os.system("mv {output}/{sampleid}/HC*.vcf {output}/HC/".format(sampleid=sampleid, output=args.outputfolder))
    os.system("mv {output}/{sampleid}/UMCU*.vcf {output}/UMCU/".format(sampleid=sampleid, output=args.outputfolder))
    os.system("rm -r {output}/{sample}".format(output=args.outputfolder, sample=sampleid))

    return [sampleid, bamfile, refset]


def parse_ped(ped_file):
    samples = {}  # 'sample_id': {'family': 'fam_id', 'parents': ['sample_id', 'sample_id']}
    for line in ped_file:
        ped_data = line.strip().split()
        family, sample, father, mother, sex, phenotype = ped_data

        # Create samples
        if sample not in samples:
            samples[sample] = {'family': family, 'parents': [], 'children': []}
        if father != '0' and father not in samples:
            samples[father] = {'family': family, 'parents': [], 'children': []}
        if mother != '0' and mother not in samples:
            samples[mother] = {'family': family, 'parents': [], 'children': []}

        # Save sample relations
        if father != '0':
            samples[sample]['parents'].append(father)
            samples[father]['children'].append(sample)
        if mother != '0':
            samples[sample]['parents'].append(mother)
            samples[mother]['children'].append(sample)
    return samples


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
        bams = (
            set(glob.glob("{}/**/*.realigned.bam".format(args.inputfolder), recursive=True))
            - set(glob.glob("{}/[eE]xome[dD]epth*/**/*.realigned.bam".format(args.inputfolder), recursive=True))
        )
    elif args.pipeline == "nf":
        bams = glob.glob("{}/bam_files/**/*.bam".format(args.inputfolder), recursive=True)
    print("Number of BAM files = "+str(len(bams)))

    if bams:
        os.system("echo \"{user} {today}\tExomeDepth reanalysis performed\" >> {inputfolder}/logbook.txt".format(
            user=user, today=today, inputfolder=args.inputfolder)
        )
    else:
        sys.exit("no bam files detected")

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
    metadata_dic = {}
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
        sampleinfo = pool.starmap(exomedepth_analysis, [[bam, args, gender_dic, suffix_dic, refset_dic] for bam in bams])

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
    for item in sampleinfo:
        action = "python {0} single_igv {1} {2} {3} {4} --bam {5}".format(
            settings.igv_xml, args.outputfolder, item[0], args.runid, item[2], item[1]
        )
        os.system(action)

    """ Make family IGV session(s)"""
    families = parse_ped(open(args.pedfile, "r"))
    for item in sampleinfo:
        if len(list(set(families[item[0]]["parents"]))) == 2:
            action = "python {0} family_igv {1} {2} {3} {4} {5}".format(
                settings.igv_xml, args.outputfolder, args.pedfile, args.runid, item[0], " ".join(bams)
            )
            os.system(action)
