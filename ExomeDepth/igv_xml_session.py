#! /usr/bin/env python
import argparse
from string import Template
import settings

def make_single_igvsession(args, sample_id, statistic, igv_extension, vcf_extension):
    template_file = Template(open(args.template_xml).read())    

    # Data files XML variables
    if args.pipeline == "iap": #For analysis based on IAP
        bam = "../{0}/mapping/{1}".format(args.sampleid, args.bam)  #BAM file
        snv_vcf = "../single_sample_vcf/{0}.filtered_variants.vcf".format(args.sampleid)  #SNV VCF file
    elif args.pipeline == "nf": #For NF pipeline.
        bam = "../bam_files/{0}".format(args.bam)  #Bam file
        snv_vcf = "../single_sample_vcf/{0}_{1}.vcf".format(args.sampleid, args.runid)  #SNV VCF file
    hc_cnv_vcf = "HC/HC_{0}_{1}_{2}_{3}".format(args.refset, args.sampleid, args.runid, vcf_extension)  # HC CNV-VCF file
    baf = "../baf/{0}.igv".format(args.sampleid)  # BAF file
    igv_hc_ratio = "igv_tracks/HC_{0}_{1}_{2}_{3}".format(args.refset, args.sampleid, args.runid, igv_extension)  #CNV igv session with ratios for HC
    igv_umcu_ratio = "igv_tracks/UMCU_{0}_{1}_{2}_{3}".format(args.refset, args.sampleid, args.runid, igv_extension)  #CNV igv session with ratios for UMCU

    # ID XML variables
    session = "{0}_{1}_{2}_igv.xml".format(sample_id, statistic, args.runid)  #Session ID
    hc_cnv_vcf_id = "CNV:HC_{0}_{1}_{2}_{3}".format(args.refset, args.sampleid, args.runid, vcf_extension)  #HC CNV -VCF track id in IGV
    snv_vcf_id = "SNV:{0}".format(snv_vcf.split("/")[-1])  #SNV track ID in IGV
    igv_hc_ratio_track = "{0}_{1}_test".format(igv_hc_ratio, statistic)  #Ratio track within HC CNV igv session
    igv_hc_ratio_track_id = "Probe_ratio:HC_{0}_{1}".format(statistic, args.sampleid)  #HC ratio track id in IGV
    igv_umcu_ratio_track = "{0}_{1}_test".format(igv_umcu_ratio, statistic)  #Ratio track within UMCU CNV igv session
    igv_umcu_ratio_track_id = "Probe_ratio:UMCU_{0}_{1}".format(statistic, args.sampleid) #UMCU ratio track id in IGV
    baf_track = "{0}_baf".format(baf)  #BAF track for BAF file
    baf_track_id = "BAF:{0}".format(sample_id)  #BAF track id in IGV
    bam_coverage = "{0}_coverage".format(bam) #BAM coverage track for BAM
    bam_id = bam.split("/")[-1]  #BAM id

    # Scale XML variables
    min_axis, mid_axis, max_axis = settings.igv_settings[statistic]  #Get axis values out of settings.py
    fontsize = args.fontsize

    # Substitue variables in IGV template
    substitute_dic = {
                     'session':session,
                     'snv_vcf':snv_vcf,
                     'igv_hc_ratio':igv_hc_ratio,
                     'igv_umcu_ratio':igv_umcu_ratio,
                     'baf':baf,
                     'bam':bam,
                     'hc_cnv_vcf':hc_cnv_vcf,
                     'snv_vcf_id':snv_vcf_id,
                     'hc_cnv_vcf_id':hc_cnv_vcf_id,
                     'igv_hc_ratio_track':igv_hc_ratio_track,
                     'igv_hc_ratio_track_id':igv_hc_ratio_track_id,
                     'mid_axis':mid_axis,
                     'max_axis':max_axis,
                     'min_axis':min_axis,
                     'baf_track':baf_track,
                     'baf_track_id':baf_track_id,
                     'igv_umcu_ratio_track':igv_umcu_ratio_track,
                     'igv_umcu_ratio_track_id':igv_umcu_ratio_track_id,
                     'bam_coverage':bam_coverage,
                     'bam_id':bam_id,
                     'sample_id':sample_id,
                     'fontsize':fontsize
                     }
    new_file = template_file.substitute(substitute_dic)
    return new_file

def make_family_igvsession(args, statistic, igv_extension, vcf_extension):
    template_file = Template(open(args.template_xml).read())
    child, father, mother = args.trioid.split(",")

    # Data files XML variables
    if args.pipeline == "iap": #For analysis based on IAP
        bam_child = "../{0}/mapping/{1}".format(child, args.bam)  #BAM file
        snv_vcf_child = "../single_sample_vcf/{0}.filtered_variants.vcf".format(child)  #SNV VCF file child
        snv_vcf_father = "../single_sample_vcf/{0}.filtered_variants.vcf".format(father)  #SNV VCF file father
        snv_vcf_mother = "../single_sample_vcf/{0}.filtered_variants.vcf".format(mother)  #SNV VCF file mother

    elif args.pipeline == "nf": #For NF pipeline.
        bam_child = "../bam_files/{0}".format(args.bam)  #Bam file
        snv_vcf_child = "../single_sample_vcf/{0}_{1}.vcf".format(child, args.runid)  #SNV VCF file child
        snv_vcf_father = "../single_sample_vcf/{0}_{1}.vcf".format(father, args.runid)  #SNV VCF file father
        snv_vcf_mother = "../single_sample_vcf/{0}_{1}.vcf".format(mother, args.runid)  #SNV VCF file mother

    hc_cnv_vcf_child = "HC/HC_{0}_{1}_{2}_{3}".format(args.refset, child, args.runid, vcf_extension)  # HC CNV-VCF file child
    hc_cnv_vcf_father = "HC/HC_{0}_{1}_{2}_{3}".format(args.refset, father, args.runid, vcf_extension)  # HC CNV-VCF file father
    hc_cnv_vcf_mother = "HC/HC_{0}_{1}_{2}_{3}".format(args.refset, mother, args.runid, vcf_extension)  # HC CNV-VCF file mother
    hc_ratio_child = "igv_tracks/HC_{0}_{1}_{2}_{3}".format(args.refset, child, args.runid, igv_extension)  #Child CNV igv session with ratios for HC
    hc_ratio_father = "igv_tracks/HC_{0}_{1}_{2}_{3}".format(args.refset, father, args.runid, igv_extension)  #Father CNV igv session with ratios for HC
    hc_ratio_mother = "igv_tracks/HC_{0}_{1}_{2}_{3}".format(args.refset, mother, args.runid, igv_extension)  #Mother CNV igv session with ratios for HC
    igv_hc_ratio = "igv_tracks/HC_{0}_{1}_{2}_{3}".format(args.refset, child, args.runid, igv_extension)  #CNV child igv session with ratios for HC
    upd_trio = "../upd/{0}".format(args.family_upd_file)  #UPD file trio
    baf_child = "../baf/{0}.igv".format(child)  # BAF file child
 
    # ID XML variables
    session = "{0}_{1}_{2}_igv.xml".format(args.familyid, args.runid, statistic)  #Session ID
    snv_vcf_child_id = "SNV_child:{0}".format(snv_vcf_child.split("/")[-1])  #SNV track child ID in IGV
    snv_vcf_father_id = "SNV_father:{0}".format(snv_vcf_father.split("/")[-1])  #SNV track father ID in IGV
    snv_vcf_mother_id = "SNV_mother:{0}".format(snv_vcf_mother.split("/")[-1])  #SNV track mother ID in IGV
    hc_cnv_vcf_child_id = "CNV_child:HC_{0}_{1}_{2}_{3}".format(args.refset, child, args.runid, vcf_extension)  #HC CNV-VCF child track id in IGV
    hc_cnv_vcf_father_id = "CNV_father:HC_{0}_{1}_{2}_{3}".format(args.refset, father, args.runid, vcf_extension)  #HC CNV-VCF father track id in IGV
    hc_cnv_vcf_mother_id = "CNV_mother:HC_{0}_{1}_{2}_{3}".format(args.refset, mother, args.runid, vcf_extension)  #HC CNV-VCF mother track id in IGV
    hc_ratio_child_track = "{0}_{1}_test".format(hc_ratio_child, statistic)  #Ratio track child within HC CNV igv session
    hc_ratio_child_track_id = "Probe_ratio_child:HC_{0}_{1}".format(statistic, child)  #HC ratio child track id in IGV
    hc_ratio_father_track = "{0}_{1}_test".format(hc_ratio_father, statistic)  #Ratio track father within HC CNV igv session
    hc_ratio_father_track_id = "Probe_ratio_father:HC_{0}_{1}".format(statistic, father)  #HC ratio father track id in IGV
    hc_ratio_mother_track = "{0}_{1}_test".format(hc_ratio_mother, statistic)  #Ratio track mother within HC CNV igv session
    hc_ratio_mother_track_id = "Probe_ratio_mother:HC_{0}_{1}".format(statistic, mother)  #HC ratio mother track id in IGV
    baf_track_child = "{0}_baf".format(baf_child)  #BAF track child for BAF file
    baf_track_child_id = "BAF_child:{0}".format(child)  #BAF track child id in IGV
    bam_coverage = "{0}_coverage".format(bam_child) #BAM coverage track child for BAM
    bam_id = bam_child.split("/")[-1]  #BAM child id

    # Scale XML variables
    min_axis, mid_axis, max_axis = settings.igv_settings[statistic]  #Get axis values out of settings.py
    fontsize = args.fontsize

    # Substitue variables in IGV template
    substitute_dic = {
                     'snv_vcf_child':snv_vcf_child,
                     'snv_vcf_father':snv_vcf_father,
                     'snv_vcf_mother':snv_vcf_mother,
                     'hc_cnv_vcf_child':hc_cnv_vcf_child,
                     'hc_cnv_vcf_father':hc_cnv_vcf_father,
                     'hc_cnv_vcf_mother':hc_cnv_vcf_mother, 
                     'hc_ratio_child':hc_ratio_child,
                     'hc_ratio_father':hc_ratio_father,
                     'hc_ratio_mother':hc_ratio_mother,
                     'upd':upd_trio,
                     'baf':baf_child,
                     'bam':bam_child,
                     'fontsize':fontsize,
                     'mid_axis':mid_axis,
                     'max_axis':max_axis,
                     'min_axis':min_axis,
                     'session':session,
                     'snv_vcf_child_id':snv_vcf_child_id,
                     'snv_vcf_father_id':snv_vcf_father_id,
                     'snv_vcf_mother_id':snv_vcf_mother_id,
                     'hc_cnv_vcf_child_id':hc_cnv_vcf_child_id,
                     'hc_cnv_vcf_father_id':hc_cnv_vcf_father_id,
                     'hc_cnv_vcf_mother_id':hc_cnv_vcf_mother_id,
                     'hc_ratio_child_track':hc_ratio_child_track,
                     'hc_ratio_father_track':hc_ratio_father_track,
                     'hc_ratio_mother_track':hc_ratio_mother_track,
                     'hc_ratio_child_track_id':hc_ratio_child_track_id,
                     'hc_ratio_father_track_id':hc_ratio_father_track_id,
                     'hc_ratio_mother_track_id':hc_ratio_mother_track_id,
                     'baf_track_child':baf_track_child,
                     'baf_track_child_id':baf_track_child_id,
                     'bam_coverage':bam_coverage,
                     'bam_id':bam_id
                     } 
    new_file = template_file.substitute(substitute_dic)
    return new_file

def make_single(args, igv_extension, vcf_extension):
    for statistic in settings.igv_settings:
        write_file = open("{0}/{1}_{2}_{3}_igv.xml".format(args.output, args.sampleid, statistic, args.runid), "w")
        write_file.write(make_single_igvsession(args, args.sampleid, statistic, igv_extension, vcf_extension))
        write_file.close()

def make_family(args, igv_extension, vcf_extension):
    for statistic in settings.igv_settings:
        write_file = open("{0}/{1}_{2}_{3}_igv.xml".format(args.output, args.familyid, args.runid, statistic), "w")
        write_file.write(make_family_igvsession(args, statistic, igv_extension, vcf_extension))
        write_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()
    
    parser_single = subparser.add_parser('single_igv', help='Make single sample IGV session')
    parser_single.add_argument('bam', help='BAM file')
    parser_single.add_argument('output', help='Output folder')
    parser_single.add_argument('sampleid', help='Sample ID')
    parser_single.add_argument('template_xml', default=settings.template_single_xml, help='Full path to template XML [default = settting.py]')
    parser_single.add_argument('refset', help='used refset (i.e. CREv2-2021-2 ')
    parser_single.add_argument('runid', help='Run ID')
    parser_single.add_argument('--pipeline', default='nf', choices=['nf', 'iap'], help='pipeline used for sample processing (nf = nexflow, IAP = illumina analysis pipeline')
    parser_single.add_argument('--vcf_filename_suffix', help='suffix that was included in the VCF filename. Do not include spaces or underscores in suffix')
    parser_single.add_argument('--fontsize', default=12, type=int, help='fontzise within IGV session for headers [default 12]')
    
    parser_single.set_defaults(func = make_single)

    parser_family = subparser.add_parser('family_igv', help='Make family IGV session')
    parser_family.add_argument('bam', help='BAM file')
    parser_family.add_argument('output', help='Output folder')
    parser_family.add_argument('familyid', help='familyid')
    parser_family.add_argument('trioid', help='comma seperated trio ID:  child, father, mother')
    parser_family.add_argument('family_upd_file', help='input family UPD file')
    parser_family.add_argument('template_xml',  default=settings.template_family_xml, help='Full path to template XML [default = settting.py]')
    parser_family.add_argument('refset', help='used refset (i.e. CREv2-2021-2 ')
    parser_family.add_argument('runid', help='Run ID')
    parser_family.add_argument('--pipeline', default='nf', choices=['nf', 'iap'], help='pipeline used for sample processing (nf = nexflow, IAP = illumina analysis pipeline')
    parser_family.add_argument('--vcf_filename_suffix', help='suffix that was included in the VCF filename. Do not include spaces or underscores in suffix')
    parser_family.add_argument('--fontsize', default=12, type=int, help='fontzise within IGV session for headers [default 12]')
    parser_family.set_defaults(func = make_family)
    args = parser.parse_args()

    igv_extension = "ref.igv"
    vcf_extension = "exome_calls.vcf"
    if args.vcf_filename_suffix:
        igv_extension = "{}_ref.igv".format(args.vcf_filename_suffix)
        vcf_extension = "exome_calls_{}.vcf".format(args.vcf_filename_suffix)
    args.func(args, igv_extension, vcf_extension)
