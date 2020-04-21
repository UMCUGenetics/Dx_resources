#! /usr/bin/env python3
from string import Template
import argparse
import settings

def make_igvsession(template_file, igv_ed_umcu, igv_ed_hc, bam, vcf_hc, sample_id, vcf_SNV, axis, statistic):
    new_session = "{0}_{1}_igv.xml".format(sample_id, statistic)
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
                    'vcf_SNV' : vcf_SNV, 'min_axis' : min_axis, 'mid_axis' : mid_axis,'max_axis' : max_axis,
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
    args = parser.parse_args()

    igv_ed_umcu = "igv_tracks/UMCU_{1}_{0}_ref.igv".format(args.bam, args.refdate)
    igv_ed_hc = "igv_tracks/HC_{1}_{0}_ref.igv".format(args.bam, args.refdate)
    bam_id = "../bam_files/{0}".format(args.bam)
    vcf_hc = "HC/HC_{0}_{1}_exome_calls.vcf".format(args.refdate, args.bam)
    vcf_SNV = "../single_sample_vcf/{0}_{1}.vcf".format(args.sampleid,args.runid)
    igv_settings = settings.igv_settings
    for statistic in igv_settings:
        write_file = open("{0}/{1}_{2}_igv.xml".format(args.output, args.sampleid, statistic), "w")
        write_file.write(make_igvsession(args.template,igv_ed_umcu, igv_ed_hc, bam_id, vcf_hc, args.sampleid, vcf_SNV, igv_settings[statistic], statistic))
        write_file.close()
