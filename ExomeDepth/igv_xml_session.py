#! /usr/bin/env python3

import argparse
import os
from string import Template

import settings


def parse_ped(ped_file):
    samples = {}  # 'sample_id': {'family': 'fam_id', 'parents': ['sample_id', 'sample_id']}
    with open(ped_file, 'r') as ped_file_lines:
        for line in ped_file_lines:
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


def parse_reanalysis_file(reanalysis_file):
    reanalysis = {}
    for line in reanalysis_file:
        sampleid, refset, tag = (
            lambda sampleid, refset, tag = None:(sampleid, refset, tag)
        )(
            *line.strip().split()
        )

        reanalysis[sampleid] = {"refset": refset, "tag": tag}
    return reanalysis


def make_single_igv_file(args, sample_id, statistic, igv_extension, vcf_extension, bam):
    template_file = Template(open(args.template_xml).read())

    """ Data files XML variables"""
    if args.pipeline == "iap":  # For analysis based on IAP
        bam_path = "../{0}/mapping/{1}".format(args.sampleid, bam)  # BAM path
        snv_vcf = "../single_sample_vcf/{0}.filtered_variants.vcf".format(args.sampleid)  # SNV VCF file
    elif args.pipeline == "nf":  # For NF pipeline.
        bam_path = "../bam_files/{0}".format(bam)  # BAM path
        snv_vcf = "../single_sample_vcf/{0}_{1}.vcf".format(args.sampleid, args.runid)  # SNV VCF file

    min_axis, mid_axis, max_axis = settings.igv_settings[statistic]  # Get axis values out of settings.py

    baf = "../baf/{0}_baf.igv".format(args.sampleid)
    igv_hc_ratio = "igv_tracks/HC_{0}_{1}_{2}_{3}".format(
        args.refset, args.sampleid, args.runid, igv_extension
    )  # CNV igv session with ratios for HC
    igv_umcu_ratio = "igv_tracks/UMCU_{0}_{1}_{2}_{3}".format(
        args.refset, args.sampleid, args.runid, igv_extension
    )  # CNV igv session with ratios for UMCU

    """ Substitue variables in IGV template"""
    substitute_dic = {
        'session': "{0}_{1}_{2}_igv.xml".format(sample_id, statistic, args.runid),
        'snv_vcf': snv_vcf,
        'igv_hc_ratio': igv_hc_ratio,
        'igv_umcu_ratio': igv_umcu_ratio,
        'baf': baf,
        'bam_path': bam_path,
        'hc_cnv_vcf': "HC/HC_{0}_{1}_{2}_{3}".format(args.refset, args.sampleid, args.runid, vcf_extension),  # HC CNV-VCF file
        'snv_vcf_id': "SNV/MNV:{0}".format(os.path.basename(snv_vcf)),  # SNV/MNV track ID in IGV
        'hc_cnv_vcf_id': "CNV:HC_{0}_{1}_{2}_{3}".format(
            args.refset, args.sampleid, args.runid, vcf_extension
        ),  # HC CNV -VCF track id in IGV
        'igv_hc_ratio_track': "{0}_{1}_test".format(igv_hc_ratio, statistic),  # Ratio track within HC CNV igv session
        'igv_hc_ratio_track_id': "Probe_ratio:HC_{0}_{1}".format(statistic, args.sampleid),  # HC ratio track id in IGV
        'mid_axis': mid_axis,
        'max_axis': max_axis,
        'min_axis': min_axis,
        'baf_track': "{0}_baf".format(baf),  # BAF track for BAF file
        'baf_track_id': "BAF:{0}".format(sample_id),  # BAF track id in IGV
        'igv_umcu_ratio_track': "{0}_{1}_test".format(igv_umcu_ratio, statistic),  # Ratio track within UMCU CNV igv session
        'igv_umcu_ratio_track_id': "Probe_ratio:UMCU_{0}_{1}".format(statistic, args.sampleid),  # UMCU ratio track id in IGV
        'bam_coverage': "{0}_coverage".format(bam_path),  # BAM coverage track for BAM
        'bam_id': bam,
        'sample_id': sample_id,
        'fontsize': args.fontsize
    }
    vcf_file = template_file.substitute(substitute_dic)
    return vcf_file


def make_single_igv_session(args):
    if args.bam:  # Use full path
        bam = args.bam
    else:  # Use relative path
        if args.pipeline == "iap":  # For analysis based on IAP
            bam = "{0}_dedup.realigned.bam".format(args.sampleid.split(",")[0])  # BAM file
        elif args.pipeline == "nf":  # For NF pipeline.
            bam = "{0}.bam".format(args.sampleid.split(",")[0])  # Bam file

    igv_extension = "ref.igv"
    vcf_extension = "exome_calls.vcf"

    if args.reanalysis:
        reanalysis = parse_reanalysis_file(args.reanalysis)
        if reanalysis[args.sampleid]['tag']:
            igv_extension = "{}_ref.igv".format(settings.reanalysis_dic[reanalysis[args.sampleid]['tag']][1])
            vcf_extension = "exome_calls_{}.vcf".format(settings.reanalysis_dic[reanalysis[args.sampleid]['tag']][1])

    for statistic in settings.igv_settings:
        write_file = open("{0}/{1}_{2}_{3}_igv.xml".format(
            args.output, args.sampleid, statistic, args.runid), "w"
        )
        write_file.write(make_single_igv_file(
            args, args.sampleid, statistic, igv_extension, vcf_extension, bam)
        )
        write_file.close()


def get_file(pattern, file_paths):
    for file_path in file_paths:
        if pattern in file_path:
            return os.path.basename(file_path)


def make_family_igv_file(args, familyid, child, father, mother, statistic):
    template_file = Template(open(args.template_xml).read())

    upd_file = get_file(child, args.upd_files)
    baf_file = get_file(child, args.baf_files)
    bam_file_child = get_file(child, args.bam_files)
    hc_ratio_child_file = get_file(child, args.igv_files)
    hc_ratio_father_file = get_file(father, args.igv_files)
    hc_ratio_mother_file = get_file(mother, args.igv_files)
    snv_vcf_child_file = get_file(child, args.snv_vcf_files)
    snv_vcf_father_file = get_file(father, args.snv_vcf_files)
    snv_vcf_mother_file = get_file(mother, args.snv_vcf_files)
    hc_cnv_vcf_child_file = get_file(child, args.cnv_vcf_files)
    hc_cnv_vcf_father_file = get_file(father, args.cnv_vcf_files)
    hc_cnv_vcf_mother_file = get_file(mother, args.cnv_vcf_files)

    #  Scale XML variables
    min_axis, mid_axis, max_axis = settings.igv_settings[statistic]  # Get axis values out of settings.py

    #  Substitue variables in IGV template
    substitute_dic = {
        "snv_vcf_child": f"../single_sample_vcf/{snv_vcf_child_file}",
        "snv_vcf_father": f"../single_sample_vcf/{snv_vcf_father_file}",
        "snv_vcf_mother": f"../single_sample_vcf/{snv_vcf_mother_file}",
        "snv_vcf_child_id": f"SNV/MNV_child:{snv_vcf_child_file}",
        "snv_vcf_father_id": f"SNV/MNV_father:{snv_vcf_father_file}",
        "snv_vcf_mother_id": f"SNV/MNV_mother:{snv_vcf_mother_file}",
        "hc_cnv_vcf_child": f"HC/{hc_cnv_vcf_child_file}",
        "hc_cnv_vcf_father": f"HC/{hc_cnv_vcf_father_file}",
        "hc_cnv_vcf_mother": f"HC/{hc_cnv_vcf_mother_file}",
        "hc_cnv_vcf_child_id": f"CNV_child:{hc_cnv_vcf_child_file}",
        "hc_cnv_vcf_father_id": f"CNV_father:{hc_cnv_vcf_father_file}",
        "hc_cnv_vcf_mother_id": f"CNV_mother:{hc_cnv_vcf_mother_file}",
        "hc_ratio_child": f"igv_tracks/{hc_ratio_child_file}",
        "hc_ratio_father": f"igv_tracks/{hc_ratio_father_file}",
        "hc_ratio_mother": f"igv_tracks/{hc_ratio_mother_file}",
        "hc_ratio_child_track": f"igv_tracks/{hc_ratio_child_file}_{statistic}_test",
        "hc_ratio_father_track": f"igv_tracks/{hc_ratio_father_file}_{statistic}_test",
        "hc_ratio_mother_track": f"igv_tracks/{hc_ratio_mother_file}_{statistic}_test",
        "hc_ratio_child_track_id": f"Probe_ratio_child:HC_{statistic}_{child}",
        "hc_ratio_father_track_id": f"Probe_ratio_father:HC_{statistic}_{father}",
        "hc_ratio_mother_track_id": f"Probe_ratio_mother:HC_{statistic}_{mother}",
        "upd": f"../upd/{upd_file}",
        "upd_track": f"../upd/{upd_file}_inheritence",
        "upd_track_id": f"Mendelian_violations_{familyid}",
        "baf": f"../baf/{baf_file}",
        "baf_track_child": f"../baf/{baf_file}_baf",
        "baf_track_child_id": f"BAF_child:{child}",
        "bam": f"../bam_files/{bam_file_child}",
        "bam_coverage": f"../bam_files/{bam_file_child}_coverage",
        "bam_id": bam_file_child,
        "fontsize": args.fontsize,
        "mid_axis": mid_axis,
        "max_axis": max_axis,
        "min_axis": min_axis,
        "session": f"{familyid}_{args.runid}_{statistic}_igv.xml"
    }
    new_file = template_file.substitute(substitute_dic)
    return new_file


def make_family_igv_session(args):
    samples = parse_ped(args.ped_file)
    for sample in samples:
        if len(samples[sample]['parents']) == 2:  # Sample = child both with parents
            familyid = samples[sample]['family']
            father = samples[sample]['parents'][0]
            mother = samples[sample]['parents'][1]
            joined_bams = "".join(args.bam_files)
            if sample in joined_bams and father in joined_bams and mother in joined_bams:  
                # Makes sure sample, father, and mother are in within analysis, not only in PED file
                for statistic in settings.igv_settings:
                    session = "FAM{0}_{1}_{2}_{3}_igv.xml".format(familyid, sample, statistic, args.runid)
                    with open(f"{args.output}/{session}", "w") as family_file:
                        family_file.write(make_family_igv_file(args, familyid, sample, father, mother, statistic))


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
        '--reanalysis',
        type=argparse.FileType('r'),
        help='Tab delimited file with SampleID, RefsetID, and optional reanalysis female/male mode'
    )
    parser_single.add_argument(
        '--fontsize', default=10, type=int,
        help='fontzise within IGV session for headers [default 12]'
    )
    parser_single.set_defaults(func=make_single_igv_session)

    parser_family = subparser.add_parser('family_igv', help='Make family IGV session')
    parser_family.add_argument('output', help='Output folder')
    parser_family.add_argument('ped_file', help='full path to ped_file')
    parser_family.add_argument('runid', help='Run ID')
    parser_family.add_argument(
        '--bam_files',
        nargs='+',
        help='full path to BAM files of all family members including child and parents (space seperated)'
    )
    parser_family.add_argument(
        '--snv_vcf_files',
        nargs='+',
        help='full path to GATK snv VCF files of all family members including child and parents (space seperated)'
    )
    parser_family.add_argument(
        '--cnv_vcf_files',
        nargs='+',
        help='full path to CNV VCF files of all family members including child and parents (space seperated)'
    )
    parser_family.add_argument(
        '--igv_files',
        nargs='+',
        help='full path to igv tracks all family members including child and parents (space seperated)'
    )
    parser_family.add_argument(
        '--upd_files',
        nargs='+',
        help='full path to UPD tracks of all family members including child and parents (space seperated)'
    )
    parser_family.add_argument(
        '--baf_files',
        nargs='+',
        help='full path to BAF tracks of all family members including child and parents (space seperated)'
    )
    parser_family.add_argument(
        '--reanalysis',
        type=argparse.FileType('r'),
        help='Tab delimited file with SampleID, RefsetID, and optional reanalysis female/male mode'
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
        '--fontsize', default=settings.fontsize, type=int,
        help='fontzise within IGV session for headers [default 12]'
    )
    parser_family.set_defaults(func=make_family_igv_session)
    args = parser.parse_args()
    args.func(args)
