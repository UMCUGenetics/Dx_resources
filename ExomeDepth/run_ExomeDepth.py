#! /usr/bin/env python3

import os
import sys
import subprocess
import re
import argparse
import glob
from multiprocessing import Pool
import pysam
import database.functions
import utils.utils
from igv_xml_session import split_fslash
import settings


def valid_read(read):
    """Check if a read is properly mapped."""
    if (read.mapping_quality >= 20 and read.reference_end and read.reference_start):
        return True
    else:
        return False


def get_gender(bam):
    """Determine chrY ratio based on read count in bam (excl PAR)."""
    with pysam.AlignmentFile(bam, "rb") as workfile:
        yreads = float(sum([valid_read(read) for read in workfile.fetch(region=settings.gender_determination_locus_y)]))
        total = float(workfile.mapped)
    yratio = float("%.2f" % ((yreads / total) * 100))
    if yratio <= float(settings.gender_determination_y_ratio[0]):
        return "female"
    elif yratio >= float(settings.gender_determination_y_ratio[1]):
        return "male"
    else:
        if re.search('[C|P]M', split_fslash(bam)):
            print("Sample {0} has a unknown gender based on chrY reads, but resolved as male based on sampleID".format(
                split_fslash(bam))
            )
            return "male"
        elif re.search('[C|P]F', split_fslash(bam)):
            print("Sample {0} has a unknown gender based on chrY reads, but resolved as female based on sampleID".format(
                split_fslash(bam))
            )
            return "female"
        else:
            sys.exit("Sample {0} has a unknown gender and will not be analysed".format(split_fslash(bam)))


def multiprocess_ref(mp_list):

    action = (
        "export SINGULARITYENV_R_LIBS={r_libs} && singularity exec -B {singularity_mnt} "
        "{singularity_container} Rscript {refscript} {folder}/ {folder}/{outputid} {targetbed} {refgenome} {exonbed}\n".format(
            r_libs=settings.r_library_path,
            singularity_mnt=settings.singularity_mount_path,
            singularity_container=settings.singularity_r_container,
            refscript=settings.create_refset_r,
            folder=mp_list[0],
            outputid=mp_list[1],
            targetbed=mp_list[2],
            refgenome=settings.reference_genome,
            exonbed=mp_list[3]
            )
        )
    os.system(action)


def make_refset(args):

    """Log all settings in setting.log file"""
    log_file = "{output}/settings.log".format(
        output=args.output
        )
    log_setting_file = open(log_file, "w")
    options = vars(args)
    for item in options:
        log_setting_file.write("{0}\t{1}\n".format(str(item), str(options[item])))
    for item in dir(settings):
        if "__" not in item:
            log_setting_file.write("{0}\t{1}\n".format(item, str(repr(eval("settings.%s" % item)))))
    log_setting_file.close()

    """Make new reference set."""
    analysis = settings.analysis
    output_folder = format(os.path.abspath(args.output))

    bams = glob.glob("{}/**/*.bam".format(args.inputfolder), recursive=True)
    print("Number of BAM files detected = {0}".format(len(bams)))

    """Get gender from chrY read count ratio."""
    ref_gender_dic = {}  # Dictionary with gender of each sample
    for bam in bams:
        gender = get_gender(bam)
        if gender != "unknown":
            if gender not in ref_gender_dic:
                ref_gender_dic[gender] = [bam]
            else:
                ref_gender_dic[gender] += [bam]
        else:
            print("Sample {0} has unknown gender and is removed from analysis".format(bam))

    """Make folder per gender + analysis, and soflink BAMs in these folders."""

    mp_list = []
    for model in analysis:
        for item in ref_gender_dic:
            folder = "{0}/{1}_{2}_{3}".format(output_folder, model, item, split_fslash(str(args.output).rstrip("/")))
            output_id = "{0}_{1}_{2}.EDref".format(model, item, args.prefix)
            os.mkdir(folder)
            for bam in ref_gender_dic[item]:
                os.symlink("{bam}* {folder}".format(bam=bam, folder=folder))
            mp_list += [[folder, output_id, analysis[model]["target_bed"], analysis[model]["exon_bed"]]]

    with Pool(processes=int(args.simjobs)) as pool:
        pool.map(multiprocess_ref, mp_list, 1)


def multiprocess_call(multiprocess_list):

    """Log all settings in setting.log file"""
    setting_log_suffix = "settings.log"
    r_log_suffix = "CNV.log"
    r_igv_suffix = "ref.igv"
    if args.vcf_filename_suffix:
        setting_log_suffix = "{0}_{1}".format(args.vcf_filename_suffix, setting_log_suffix)
        r_log_suffix = "{0}_{1}".format(args.vcf_filename_suffix, r_log_suffix)
        r_igv_suffix = "{0}_{1}".format(args.vcf_filename_suffix, r_igv_suffix)

    log_file = "{output}/{model}_{refset}_{bam}_{run}_{setting_log_suffix}".format(
        output=multiprocess_list["output_folder"],
        model=multiprocess_list["model"],
        refset=multiprocess_list["refset"],
        bam=args.sample,
        run=args.run,
        setting_log_suffix=setting_log_suffix
        )
    log_setting_file = open(log_file, "w")
    options = vars(args)
    for item in options:
        log_setting_file.write("{0}\t{1}\n".format(str(item), str(options[item])))
    for item in dir(settings):
        if "__" not in item:
            log_setting_file.write("{0}\t{1}\n".format(item, str(repr(eval("settings.%s" % item)))))
    log_setting_file.close()

    """Perform ExomeDepth analysis"""
    refset_R = "{refset_dir}/{model}_{gender}_{refset}.EDref".format(
        refset_dir=settings.refset_dir,
        model=multiprocess_list["model"],
        gender=multiprocess_list["gender"],
        refset=multiprocess_list["refset"]
        )

    action = (
        "export SINGULARITYENV_R_LIBS={r_libs} && "
        "singularity exec -B {singularity_mnt} {singularity_container} "
        "Rscript {ed_r} {refset_R} {target_bed} {refgenome} {exon_bed} "
        "{prob} {bam} {model} {refset} {expected} {run} {r_log_suffix} "
        "{r_igv_suffix} {sampleid}").format(
            r_libs=settings.r_library_path,
            singularity_mnt=settings.singularity_mount_path,
            singularity_container=settings.singularity_r_container,
            ed_r=settings.call_cnv_r,
            refset_R=refset_R,
            target_bed=multiprocess_list["target_bed"],
            refgenome=settings.reference_genome,
            exon_bed=multiprocess_list["exon_bed"],
            prob=settings.probability[str(multiprocess_list["model"])],
            bam=multiprocess_list["bam"],
            model=multiprocess_list["model"],
            refset=multiprocess_list["refset"],
            expected=args.expectedCNVlength,
            run=args.run,
            r_log_suffix=r_log_suffix,
            r_igv_suffix=r_igv_suffix,
            sampleid=args.sample
        )
    os.system(action)

    """Perform csv to vcf conversion """
    inputcsv = "{0}/{1}_{2}_{3}_{4}_exome_calls.csv".format(
        multiprocess_list["output_folder"],
        multiprocess_list["model"],
        multiprocess_list["refset"],
        args.sample,
        args.run
        )

    action = (
        "python {csv2vcf} {inputcsv} {refset} {calling_model} {gender} "
        " {sampleid} {template} {runid} ").format(
            csv2vcf=settings.csv2vcf,
            inputcsv=inputcsv,
            refset=multiprocess_list["refset"],
            calling_model=multiprocess_list["calling_model"],
            gender=multiprocess_list["gender"],
            sampleid=args.sample,
            template=settings.vcf_template,
            runid=args.run
            )

    if args.vcf_filename_suffix:
        action = "{action} --vcf_filename_suffix {vcf_filename_suffix}".format(
            action=action, vcf_filename_suffix=args.vcf_filename_suffix
        )

    os.system(action)


def call_cnv(args):

    """Call CNV from BAMs"""
    bam = format(os.path.abspath(args.inputbam))
    output_folder = format(os.path.abspath(args.output))
    analysis = settings.analysis

    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    """Determine gender"""
    if args.refset_gender:  # Used gender if used as input parameter.
        gender = args.refset_gender
    else:  # Otherwise determine based on chrY count
        gender = get_gender(bam)

    if args.refset:  # Do not query database to query refset
        refset = args.refset
    else:
        """ Add sample to database if not present, or query refset from db if present """
        refset = database.functions.add_sample_to_db_and_return_refset_bam(bam, settings.refset, False)

    multiprocess_list = []
    for model in analysis:
        multiprocess_list.append({
            "model": model,
            "output_folder": output_folder,
            "gender": gender,
            "target_bed": analysis[model]["target_bed"],
            "exon_bed": analysis[model]["exon_bed"],
            "calling_model": analysis[model]["calling_model"],
            "bam": bam,
            "refset": refset
        })

    with Pool(processes=int(args.simjobs)) as pool:
        pool.map(multiprocess_call, multiprocess_list, 1)

    """Make log for stats of each model """
    for model in analysis:

        stats_log_suffix = "_"
        if args.vcf_filename_suffix:
            stats_log_suffix = "{0}{1}_".format(stats_log_suffix, args.vcf_filename_suffix)

        sample_model_log = open("{output}/{model}_{sample}{stats_log_suffix}stats.log".format(
            output=args.output,
            model=model,
            stats_log_suffix=stats_log_suffix,
            sample=args.sample
            ), "w")

        """ Get stats from VCF """
        vcf_suffix = "exome_calls"
        if args.vcf_filename_suffix:
            vcf_suffix = "{}_{}".format(vcf_suffix, args.vcf_filename_suffix)

        vcf = "{output}/{model}_{refset}_{sample}_{run}_{vcf_suffix}.vcf".format(
            output=args.output,
            model=model,
            refset=refset,
            sample=args.sample,
            run=args.run,
            vcf_suffix=vcf_suffix
            )

        stats = (subprocess.getoutput("tail -n1 {}".format(vcf)).split()[-1]).split(":")
        correlation, del_dup_ratio, number_calls = float(stats[4]), float(stats[8]), int(stats[9])

        qc_status = ""
        if args.qc_stats:
            if(correlation < float(settings.correlation) or
               number_calls < int(settings.number_calls[0]) or
               number_calls > int(settings.number_calls[1]) or
               del_dup_ratio < float(settings.del_dup_ratio[0]) or
               del_dup_ratio > float(settings.del_dup_ratio[1])):
                qc_status = "{qc_status}\tWARNING:QC_FAIL".format(qc_status=qc_status)
        if args.vcf_filename_suffix:
            qc_status = "{qc_status}\tWARNING:{qc_suffix}".format(qc_status=qc_status, qc_suffix=args.vcf_filename_suffix)

        sample_model_log.write(
            "{sample}\t{model}\t{refset}\t{correlation}\t{del_dup_ratio}\t{number_calls}{qc_status}\n".format(
                sample=args.sample,
                model=model,
                refset=refset,
                correlation=correlation,
                del_dup_ratio=del_dup_ratio,
                number_calls=number_calls,
                qc_status=qc_status
            )
        )

        sample_model_log.close()


def call_exomedepth_summary(args):
    utils.utils.exomedepth_summary(args.exomedepth_logs, args.print_stdout)


def call_detect_merge(args):
    utils.utils.detect_merge(args.inputfolder, args.outputfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()

    parser_refset = subparser.add_parser('makeref', help='Make ExomeDepth reference set')
    parser_refset.add_argument('output', help='Output folder for reference set files')
    parser_refset.add_argument('inputfolder', help='Input folder containing BAM files')
    parser_refset.add_argument('prefix', help='Prefix for reference set (e.g. Jan2020)')
    parser_refset.add_argument(
        '--simjobs', default=4,
        help='number of simultaneous samples to proces. Note: make sure similar threads are reseved in session! [default = 4]'
    )
    parser_refset.set_defaults(func=make_refset)

    parser_cnv = subparser.add_parser('callcnv', help='Call CNV with ExomeDepth basedon BAM file')
    parser_cnv.add_argument('output', help='output folder for CNV calling')
    parser_cnv.add_argument('inputbam', help='Input BAM file')
    parser_cnv.add_argument('run', help='Name of the run')
    parser_cnv.add_argument('sample', help='Sample name')
    parser_cnv.add_argument('--refset', help='Reference set to be used (e.g. CREv2-2021-2).')
    parser_cnv.add_argument(
        '--template', default=settings.template_single_xml,
        help='Template XML for single sample IGV session. Default = template_single_xml in settings.py'
    )
    parser_cnv.add_argument(
        '--pipeline', default='nf', choices=['nf', 'iap'],
        help='pipeline used for sample processing (nf = nexflow (default), IAP = illumina analysis pipeline'
    )
    parser_cnv.add_argument(
        '--simjobs', default=2,
        help='number of simultaneous samples to proces. Note: make sure similar threads are reseved in session! [default = 2]'
    )
    parser_cnv.add_argument(
        '--refset_gender', choices=['male', 'female'],
        help='force specific use of female/male reference set in analysis'
    )
    parser_cnv.add_argument(
        '--vcf_filename_suffix',
        help='suffix to be included in VCF filename. Do not include spaces or underscores in suffix'
    )
    parser_cnv.add_argument(
        '--expectedCNVlength', default=settings.expectedCNVlength,
        help='expected CNV length (basepairs) taken into account by ExomeDepth [default expectedCNVlength in settings.py]'
    )
    parser_cnv.add_argument(
        '--qc_stats', action='store_true',
        help='switch on QC check for exomedepth VCF stats (default = off)'
    )
    parser_cnv.set_defaults(func=call_cnv)

    parser_summary = subparser.add_parser('summary', help='Make ExomeDepth summary file')
    parser_summary.add_argument('exomedepth_logs', type=argparse.FileType('r'), nargs='*', help='Exomedepth log files')
    parser_summary.add_argument('--print_stdout', default=True, help='print output in stdout [default = True]')
    parser_summary.set_defaults(func=call_exomedepth_summary)

    parser_merge = subparser.add_parser('identify_merge', help='Identify merge samples')
    parser_merge.add_argument('inputfolder', help='input folder which included BAM files')
    parser_merge.add_argument('outputfile', help='output filename of identified merge samples')
    parser_merge.set_defaults(func=call_detect_merge)

    args = parser.parse_args()
    args.func(args)
