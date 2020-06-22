#! /usr/bin/env python
import argparse
from string import Template
import settings

def make_igvsession(igv_ed_umcu, igv_ed_hc, bam, vcf_hc, sample_id, vcf_SNV, axis, statistic):
    template_file = Template(open(settings.template_xml).read())
    new_session = "{0}_{1}_{2}_igv.xml".format(sample_id, statistic, args.runid)
    igv_ed_hc_test = "{0}_{1}_test".format(igv_ed_hc, statistic)
    igv_ed_umcu_test = "{0}_{1}_test".format(igv_ed_umcu, statistic)
    bam_coverage = "{0}_coverage".format(bam)
    bam_junctions = "{0}_junctions".format(bam)
    min_axis, mid_axis, max_axis = axis
    ratioid_UMCU = "{0}_UMCU".format(statistic)
    ratioid_HC = "{0}_HC".format(statistic)
    substitute_dic={'session_var' : new_session, 'igv_ed_umcu' : igv_ed_umcu, 'igv_ed_hc' : igv_ed_hc, 
                    'bam':bam, 'vcf_hc' : vcf_hc, 'sample_id' : sample_id, 'igv_ed_hc_test' : igv_ed_hc_test,
                    'igv_ed_umcu_test' : igv_ed_umcu_test, 'bam_coverage' : bam_coverage, 'bam_junctions' : bam_junctions, 
                    'vcf_SNV' : vcf_SNV, 'min_axis' : min_axis, 'mid_axis' : mid_axis, 'max_axis' : max_axis,
                    'ratioid_UMCU' : ratioid_UMCU, 'ratioid_HC': ratioid_HC 
                   }
    new_file = template_file.substitute(substitute_dic)
    return new_file

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help='BAM file')
    parser.add_argument('output', help='Output folder')
    parser.add_argument('sampleid', help='Sample ID')
    parser.add_argument('template', help='Full path to template XML')
    parser.add_argument('refdate', help='Date of the used reference set')
    parser.add_argument('runid', help='Run ID')
    parser.add_argument('--pipeline', default='nf', choices=['nf', 'iap'], help='pipeline used for sample processing (nf = nexflow, IAP = illumina analysis pipeline')
    args = parser.parse_args()

    igv_settings = settings.igv_settings
    igv_ed_umcu = "igv_tracks/UMCU_{0}_{1}_{2}_ref.igv".format(args.refdate, args.bam, args.runid)
    igv_ed_hc = "igv_tracks/HC_{0}_{1}_{2}_ref.igv".format(args.refdate, args.bam, args.runid)
    vcf_hc = "HC/HC_{0}_{1}_{2}_exome_calls.vcf".format(args.refdate, args.bam, args.runid)
    if args.pipeline == "iap": #For analysis based on IAP
        bam_id = "../{0}/mapping/{1}".format(args.sampleid, args.bam)
        vcf_SNV = "../single_sample_vcf/{0}.filtered_variants.vcf".format(args.sampleid)
    elif args.pipeline == "nf": #For NF pipeline.
        bam_id = "../bam_files/{0}.bam".format(args.sampleid)
        vcf_SNV = "../single_sample_vcf/{0}_{1}.vcf".format(args.sampleid, args.runid)
 
    for statistic in igv_settings:
        write_file = open("{0}/{1}_{2}_{3}_igv.xml".format(args.output, args.sampleid, statistic, args.runid), "w")
        write_file.write(make_igvsession(igv_ed_umcu, igv_ed_hc, bam_id, vcf_hc, args.sampleid, vcf_SNV, igv_settings[statistic], statistic))
        write_file.close()
