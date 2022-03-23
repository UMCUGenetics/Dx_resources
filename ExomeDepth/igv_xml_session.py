#! /usr/bin/env python3

import argparse
from string import Template
from exomedepth_db import add_sample_to_db_and_return_refset_bam
import settings


def parse_ped(ped_file, sampleid):
    sample_dic = {}
    for line in ped_file:
        splitline = line.split()
        if splitline[1] not in sample_dic:
            sample_dic[splitline[1]] = [splitline[0], splitline[2], splitline[3]]
    return sample_dic[sampleid][0], sampleid, sample_dic[sampleid][1], sample_dic[sampleid][2]


def make_single_igvsession(args, sample_id, statistic, igv_extension, vcf_extension, bam):
    template_file = Template(open(args.template_xml).read())

    """ Data files XML variables"""
    if args.pipeline == "iap":  # For analysis based on IAP
        bam_path = "../{0}/mapping/{1}".format(args.sampleid, bam)  # BAM path
        snv_vcf = "../single_sample_vcf/{0}.filtered_variants.vcf".format(args.sampleid)  # SNV VCF file
    elif args.pipeline == "nf":  # For NF pipeline.
        bam_path = "../bam_files/{0}".format(bam)  # BAM path
        snv_vcf = "../single_sample_vcf/{0}_{1}.vcf".format(args.sampleid, args.runid)  # SNV VCF file

    hc_cnv_vcf = "HC/HC_{0}_{1}_{2}_{3}".format(args.refset, args.sampleid, args.runid, vcf_extension)  # HC CNV-VCF file
    baf = "../baf/{0}_baf.igv".format(args.sampleid)  # BAF file
    igv_hc_ratio = "igv_tracks/HC_{0}_{1}_{2}_{3}".format(
        args.refset, args.sampleid, args.runid, igv_extension
    )  # CNV igv session with ratios for HC
    igv_umcu_ratio = "igv_tracks/UMCU_{0}_{1}_{2}_{3}".format(
        args.refset, args.sampleid, args.runid, igv_extension
    )  # CNV igv session with ratios for UMCU

    """ ID XML variables"""
    session = "{0}_{1}_{2}_igv.xml".format(sample_id, statistic, args.runid)  # Session ID
    hc_cnv_vcf_id = "CNV:HC_{0}_{1}_{2}_{3}".format(
        args.refset, args.sampleid, args.runid, vcf_extension
    )  # HC CNV -VCF track id in IGV
    snv_vcf_id = "SNV/MNV:{0}".format(snv_vcf.split("/")[-1])  # SNV/MNV track ID in IGV
    igv_hc_ratio_track = "{0}_{1}_test".format(igv_hc_ratio, statistic)  # Ratio track within HC CNV igv session
    igv_hc_ratio_track_id = "Probe_ratio:HC_{0}_{1}".format(statistic, args.sampleid)  # HC ratio track id in IGV
    igv_umcu_ratio_track = "{0}_{1}_test".format(igv_umcu_ratio, statistic)  # Ratio track within UMCU CNV igv session
    igv_umcu_ratio_track_id = "Probe_ratio:UMCU_{0}_{1}".format(statistic, args.sampleid)  # UMCU ratio track id in IGV
    baf_track = "{0}_baf".format(baf)  # BAF track for BAF file
    baf_track_id = "BAF:{0}".format(sample_id)  # BAF track id in IGV
    bam_coverage = "{0}_coverage".format(bam_path)  # BAM coverage track for BAM

    """ Scale XML variables"""
    min_axis, mid_axis, max_axis = settings.igv_settings[statistic]  # Get axis values out of settings.py
    fontsize = args.fontsize

    """ Substitue variables in IGV template"""
    substitute_dic = {
        'session': session,
        'snv_vcf': snv_vcf,
        'igv_hc_ratio': igv_hc_ratio,
        'igv_umcu_ratio': igv_umcu_ratio,
        'baf': baf,
        'bam_path': bam_path,
        'hc_cnv_vcf': hc_cnv_vcf,
        'snv_vcf_id': snv_vcf_id,
        'hc_cnv_vcf_id': hc_cnv_vcf_id,
        'igv_hc_ratio_track': igv_hc_ratio_track,
        'igv_hc_ratio_track_id': igv_hc_ratio_track_id,
        'mid_axis': mid_axis,
        'max_axis': max_axis,
        'min_axis': min_axis,
        'baf_track': baf_track,
        'baf_track_id': baf_track_id,
        'igv_umcu_ratio_track': igv_umcu_ratio_track,
        'igv_umcu_ratio_track_id': igv_umcu_ratio_track_id,
        'bam_coverage': bam_coverage,
        'bam_id': bam,
        'sample_id': sample_id,
        'fontsize': fontsize
    }
    new_file = template_file.substitute(substitute_dic)
    return new_file


def make_family_igvsession(args, statistic, igv_extension, vcf_extension, familyid, child, father, mother, child_bam, refsets):
    template_file = Template(open(args.template_xml).read())
    child_ref, father_ref, mother_ref = refsets

    """ Data files XML variables """
    if args.pipeline == "iap":  # For analysis based on IAP
        bam_child_path = "../{0}/mapping/{1}".format(child, child_bam)  # BAM path
        snv_vcf_child = "../single_sample_vcf/{0}.filtered_variuants.vcf".format(child)  # SNV VCF file child
        snv_vcf_father = "../single_sample_vcf/{0}.filtered_variants.vcf".format(father)  # SNV VCF file father
        snv_vcf_mother = "../single_sample_vcf/{0}.filtered_variants.vcf".format(mother)  # SNV VCF file mother

    elif args.pipeline == "nf":  # For NF pipeline.
        bam_child_path = "../bam_files/{0}".format(child_bam)  # BAM path
        snv_vcf_child = "../single_sample_vcf/{0}_{1}.vcf".format(child, args.runid)  # SNV VCF file child
        snv_vcf_father = "../single_sample_vcf/{0}_{1}.vcf".format(father, args.runid)  # SNV VCF file father
        snv_vcf_mother = "../single_sample_vcf/{0}_{1}.vcf".format(mother, args.runid)  # SNV VCF file mother

    hc_cnv_vcf_child = "HC/HC_{0}_{1}_{2}_{3}".format(
        child_ref, child, args.runid, vcf_extension
    )  # HC CNV-VCF file child
    hc_cnv_vcf_father = "HC/HC_{0}_{1}_{2}_{3}".format(
        father_ref, father, args.runid, vcf_extension
    )  # HC CNV-VCF file father
    hc_cnv_vcf_mother = "HC/HC_{0}_{1}_{2}_{3}".format(
        mother_ref, mother, args.runid, vcf_extension
    )  # HC CNV-VCF file mother
    hc_ratio_child = "igv_tracks/HC_{0}_{1}_{2}_{3}".format(
        child_ref, child, args.runid, igv_extension
    )  # Child CNV igv session with ratios for HC
    hc_ratio_father = "igv_tracks/HC_{0}_{1}_{2}_{3}".format(
        father_ref, father, args.runid, igv_extension
    )  # Father CNV igv session with ratios for HC
    hc_ratio_mother = "igv_tracks/HC_{0}_{1}_{2}_{3}".format(
        mother_ref, mother, args.runid, igv_extension
    )  # Mother CNV igv session with ratios for HC
    upd = "../upd/{0}_{1}.igv".format(args.runid, familyid)  # UPD file
    upd_track = "{0}_inheritence".format(upd)
    baf_child = "../baf/{0}_baf.igv".format(child)  # BAF file child

    """ ID XML variables"""
    session = "{0}_{1}_{2}_igv.xml".format(
        familyid, args.runid, statistic
    )  # Session ID
    snv_vcf_child_id = "SNV/MNV_child:{0}".format(
        snv_vcf_child.split("/")[-1]
    )  # SNV/MNV track child ID in IGV
    snv_vcf_father_id = "SNV/MNV_father:{0}".format(
        snv_vcf_father.split("/")[-1]
    )  # SNV/MNV track father ID in IGV
    snv_vcf_mother_id = "SNV/MNV_mother:{0}".format(
        snv_vcf_mother.split("/")[-1]
    )  # SNV/MNV track mother ID in IGV
    hc_cnv_vcf_child_id = "CNV_child:HC_{0}_{1}_{2}_{3}".format(
        child_ref, child, args.runid, vcf_extension
    )  # HC CNV-VCF child track id in IGV
    hc_cnv_vcf_father_id = "CNV_father:HC_{0}_{1}_{2}_{3}".format(
        father_ref, father, args.runid, vcf_extension
    )  # HC CNV-VCF father track id in IGV
    hc_cnv_vcf_mother_id = "CNV_mother:HC_{0}_{1}_{2}_{3}".format(
        mother_ref, mother, args.runid, vcf_extension
    )  # HC CNV-VCF mother track id in IGV
    hc_ratio_child_track = "{0}_{1}_test".format(
        hc_ratio_child, statistic
    )  # Ratio track child within HC CNV igv session
    hc_ratio_child_track_id = "Probe_ratio_child:HC_{0}_{1}".format(
        statistic, child
    )  # HC ratio child track id in IGV
    hc_ratio_father_track = "{0}_{1}_test".format(
        hc_ratio_father, statistic
    )  # Ratio track father within HC CNV igv session
    hc_ratio_father_track_id = "Probe_ratio_father:HC_{0}_{1}".format(
        statistic, father
    )  # HC ratio father track id in IGV
    hc_ratio_mother_track = "{0}_{1}_test".format(
        hc_ratio_mother, statistic
    )  # Ratio track mother within HC CNV igv session
    hc_ratio_mother_track_id = "Probe_ratio_mother:HC_{0}_{1}".format(
        statistic, mother
    )  # HC ratio mother track id in IGV
    baf_track_child = "{0}_baf".format(baf_child)  # BAF track child for BAF file
    baf_track_child_id = "BAF_child:{0}".format(child)  # BAF track child id in IGV
    bam_coverage = "{0}_coverage".format(bam_child_path)  # BAM coverage track child for BAM
    upd_track_id = "Mendelian_violations_{0}".format(familyid)

    """ Scale XML variables"""
    min_axis, mid_axis, max_axis = settings.igv_settings[statistic]  # Get axis values out of settings.py
    fontsize = args.fontsize

    """ Substitue variables in IGV template"""
    substitute_dic = {
        'snv_vcf_child': snv_vcf_child,
        'snv_vcf_father': snv_vcf_father,
        'snv_vcf_mother': snv_vcf_mother,
        'hc_cnv_vcf_child': hc_cnv_vcf_child,
        'hc_cnv_vcf_father': hc_cnv_vcf_father,
        'hc_cnv_vcf_mother': hc_cnv_vcf_mother,
        'hc_ratio_child': hc_ratio_child,
        'hc_ratio_father': hc_ratio_father,
        'hc_ratio_mother': hc_ratio_mother,
        'upd': upd,
        'upd_track': upd_track,
        'upd_track_id': upd_track_id,
        'baf': baf_child,
        'bam': bam_child_path,
        'fontsize': fontsize,
        'mid_axis': mid_axis,
        'max_axis': max_axis,
        'min_axis': min_axis,
        'session': session,
        'snv_vcf_child_id': snv_vcf_child_id,
        'snv_vcf_father_id': snv_vcf_father_id,
        'snv_vcf_mother_id': snv_vcf_mother_id,
        'hc_cnv_vcf_child_id': hc_cnv_vcf_child_id,
        'hc_cnv_vcf_father_id': hc_cnv_vcf_father_id,
        'hc_cnv_vcf_mother_id': hc_cnv_vcf_mother_id,
        'hc_ratio_child_track': hc_ratio_child_track,
        'hc_ratio_father_track': hc_ratio_father_track,
        'hc_ratio_mother_track': hc_ratio_mother_track,
        'hc_ratio_child_track_id': hc_ratio_child_track_id,
        'hc_ratio_father_track_id': hc_ratio_father_track_id,
        'hc_ratio_mother_track_id': hc_ratio_mother_track_id,
        'baf_track_child': baf_track_child,
        'baf_track_child_id': baf_track_child_id,
        'bam_coverage': bam_coverage,
        'bam_id': child_bam
    }
    new_file = template_file.substitute(substitute_dic)
    return new_file


def make_single(args, igv_extension, vcf_extension):
    if args.bam:  # Use full path
        bam = args.bam
    else:  # Use relative path
        if args.pipeline == "iap":  # For analysis based on IAP
            bam = "{0}_dedup.realigned.bam".format(args.sampleid.split(",")[0])  # BAM file
        elif args.pipeline == "nf":  # For NF pipeline.
            bam = "{0}.bam".format(args.sampleid.split(",")[0])  # Bam file

    for statistic in settings.igv_settings:
        write_file = open("{0}/{1}_{2}_{3}_igv.xml".format(
            args.output, args.sampleid, statistic, args.runid), "w"
        )
        write_file.write(make_single_igvsession(
            args, args.sampleid, statistic, igv_extension, vcf_extension, bam)
        )
        write_file.close()


def get_refset(bam_file, sample):
    class refset_arguments:
        bam = bam_file
        sample_id = sample
        print_refset = False
    return add_sample_to_db_and_return_refset_bam(refset_arguments)


def get_bam(sample, bam_files):
    for bam in bam_files:
        if sample in bam:
            return bam


def make_family(args, igv_extension, vcf_extension):

    familyid, child, father, mother = parse_ped(args.ped_file, args.sampleid)

    child_bam = get_bam(child, args.bam_files)
    father_bam = get_bam(father, args.bam_files)
    mother_bam = get_bam(mother, args.bam_files)

    refsets = [
        get_refset(child_bam, child),
        get_refset(father_bam, father),
        get_refset(mother_bam, mother)
    ]

    for statistic in settings.igv_settings:
        write_file = open("{0}/FAM{1}_{2}_{3}_{4}_igv.xml".format(
            args.output, familyid, child, args.runid, statistic), "w"
        )
        write_file.write(make_family_igvsession(
            args, statistic, igv_extension, vcf_extension,
            familyid, child, father, mother, child_bam, refsets
            )
        )
        write_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()

    parser_single = subparser.add_parser('single_igv', help='Make single sample IGV session')
    parser_single.add_argument('output', help='Output folder')
    parser_single.add_argument('sampleid', help='Sample ID')
    parser_single.add_argument('runid', help='Run ID')
    parser_single.add_argument('refset', help='refset of sample')
    parser_single.add_argument('--bam', help="BAM file ID")
    parser_single.add_argument(
        '--template_xml', default=settings.template_single_xml, help='Full path to template XML [default = settting.py]'
    )
    parser_single.add_argument(
        '--pipeline', default='nf', choices=['nf', 'iap'],
        help='pipeline used for sample processing (nf = nexflow, IAP = illumina analysis pipeline'
    )
    parser_single.add_argument(
        '--vcf_filename_suffix',
        help='suffix that was included in the VCF filename. Do not include spaces or underscores in suffix'
    )
    parser_single.add_argument(
        '--fontsize', default=10, type=int,
        help='fontzise within IGV session for headers [default 12]'
    )
    parser_single.set_defaults(func=make_single)

    parser_family = subparser.add_parser('family_igv', help='Make family IGV session')
    parser_family.add_argument('output', help='Output folder')
    parser_family.add_argument('ped_file', type=argparse.FileType('r'), help='full path to ped_file')
    parser_family.add_argument('runid', help='Run ID')
    parser_family.add_argument('sampleid', help='Sample ID (child)')
    parser_family.add_argument(
        'bam_files',
        nargs='+',
        help='full path to BAM files of all family members of child (space seperated)'
    )
    parser_family.add_argument(
        '--template_xml',  default=settings.template_family_xml,
        help='Full path to template XML [default = settting.py]'
    )
    parser_family.add_argument(
        '--pipeline', default='nf', choices=['nf', 'iap'],
        help='pipeline used for sample processing (nf = nexflow, IAP = illumina analysis pipeline'
    )
    parser_family.add_argument(
        '--vcf_filename_suffix',
        help='suffix that was included in the VCF filename. Do not include spaces or underscores in suffix'
    )
    parser_family.add_argument(
        '--fontsize', default=10, type=int,
        help='fontzise within IGV session for headers [default 12]'
    )
    parser_family.set_defaults(func=make_family)
    args = parser.parse_args()

    igv_extension = "ref.igv"
    vcf_extension = "exome_calls.vcf"
    if args.vcf_filename_suffix:
        igv_extension = "{}_ref.igv".format(args.vcf_filename_suffix)
        vcf_extension = "exome_calls_{}.vcf".format(args.vcf_filename_suffix)

    args.func(args, igv_extension, vcf_extension)
