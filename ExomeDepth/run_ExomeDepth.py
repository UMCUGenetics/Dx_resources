#! /usr/bin/env python3

import os
import sys
import subprocess
import re
import argparse
import glob
from multiprocessing import Pool
import pysam
from database import connect_database
from models import Sample
from exomedepth_db import create_sample
from exomedepth_db import store_sample
import settings


def valid_read(read):
    """Check if a read is properly mapped."""
    if (read.mapping_quality >= 20 and read.reference_end and read.reference_start):
        return True
    else:
        return False


def get_gender(bam):
    """Determine chrY ratio based on read count in bam (excl PAR)."""
    workfile = pysam.AlignmentFile(bam, "rb")
    yreads = float(sum([valid_read(read) for read in workfile.fetch(region=settings.gender_determination_locus_y)]))
    total = float(workfile.mapped)
    yratio = float("%.2f" % ((yreads / total) * 100))
    if yratio <= float(settings.gender_determination_y_ratio[0]):
        return "female"
    elif yratio >= float(settings.gender_determination_y_ratio[1]):
        return "male"
    else:
        if re.search('[C|P]M', bam.split("/")[-1]):
            print("Sample {0} has a unknown gender based on chrY reads, but resolved as male based on sampleID".format(
                bam.split("/")[-1])
            )
            return "male"
        elif re.search('[C|P]F', bam.split("/")[-1]):
            print("Sample {0} has a unknown gender based on chrY reads, but resolved as female based on sampleID".format(
                bam.split("/")[-1])
            )
            return "female"
        else:
            sys.exit("Sample {0} has a unknown gender and will not be analysed".format(bam.split("/")[-1]))


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
            folder = "{0}/{1}_{2}_{3}".format(output_folder, model, item, str(args.output).rstrip("/").split("/")[-1])
            output_id = "{0}_{1}_{2}.EDref".format(model, item, args.prefix)
            action = "mkdir -p {0}".format(folder)
            os.system(action)
            for bam in ref_gender_dic[item]:
                os.system("ln -sd " + str(bam) + "* " + str(folder))
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
        output=multiprocess_list[1],
        model=multiprocess_list[0],
        refset=multiprocess_list[7],
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
        model=multiprocess_list[0],
        gender=multiprocess_list[2],
        refset=multiprocess_list[7]
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
            target_bed=multiprocess_list[3],
            refgenome=settings.reference_genome,
            exon_bed=multiprocess_list[4],
            prob=settings.probability[str(multiprocess_list[0])],
            bam=multiprocess_list[5],
            model=multiprocess_list[0],
            refset=multiprocess_list[7],
            expected=args.expectedCNVlength,
            run=args.run,
            r_log_suffix=r_log_suffix,
            r_igv_suffix=r_igv_suffix,
            sampleid=args.sample
        )
    os.system(action)

    """Perform csv to vcf conversion """
    inputcsv = "{0}/{1}_{2}_{3}_{4}_exome_calls.csv".format(
        multiprocess_list[6],
        multiprocess_list[0],
        multiprocess_list[7],
        args.sample,
        args.run
        )

    action = (
        "python {csv2vcf} {inputcsv} {refset} {calling_model} {gender} "
        " {sampleid} {template} {runid} ").format(
            csv2vcf=settings.csv2vcf,
            inputcsv=inputcsv,
            refset=args.refset,
            calling_model=multiprocess_list[7],
            gender=multiprocess_list[2],
            sampleid=args.sample,
            template=settings.vcf_template,
            runid=args.run
            )

    if args.vcf_filename_suffix:
        action = "{action} --vcf_filename_suffix {vcf_filename_suffix}".format(
            action=action, vcf_filename_suffix=args.vcf_filename_suffix
        )

    os.system(action)


def add_sample_to_db(sample_id, flowcell_id, refset):
    Session = connect_database()
    with Session() as session:
        if not session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).all():
            entry = create_sample(sample_id, flowcell_id, refset)
            store_sample(Session, entry)


def query_refset(sample_id, flowcell_id):
    Session = connect_database()
    with Session() as session:
        return session.query(Sample).filter(Sample.sample == sample_id).filter(Sample.flowcell == flowcell_id).one().refset


def get_flowcelid(bam):
    workfile = pysam.AlignmentFile(bam, "rb")
    readgroups = []
    for readgroup in workfile.header['RG']:
        if readgroup['PU'] not in readgroup:
            readgroups.append(readgroup['PU'])
    readgroups = list(set(readgroups))
    flowcell_id = "_".join(readgroups)
    return flowcell_id


def call_cnv(args):

    """Call CNV from BAMs"""
    bam = format(os.path.abspath(args.inputbam))
    output_folder = format(os.path.abspath(args.output))
    analysis = settings.analysis
    if not os.path.isdir(output_folder):
        os.system("mkdir -p " + str(output_folder))

    """Determine gender"""
    if args.refset_gender:  # Used gender if used as input parameter.
        gender = args.refset_gender
    else:  # Otherwise determine based on chrY count
        gender = get_gender(bam)

    flowcell_id = get_flowcelid(bam)

    """Determine refset"""
    refset = settings.refset  # Add sample with default refset in settings.py
    """ Add sample to database if not present, and use current production refset if not present in db  """
    add_sample_to_db(args.sample, flowcell_id, refset)
    if args.refset:  # Do not use refset in database, but from argument (overrule)
        refset = args.refset
    else:  # Use refset in db
        refset = query_refset(args.sample, flowcell_id)
        # IS THIS NEEDED??? ###############################################################################

    multiprocess_list = []
    for model in analysis:
        multiprocess_list += [[(
            model, output_folder, gender, analysis[model]["target_bed"],
            analysis[model]["exon_bed"], bam, output_folder, refset
            )]
        ]

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

    args = parser.parse_args()
    args.func(args)
