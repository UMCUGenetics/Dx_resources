#! /usr/bin/env python
from string import Template
from optparse import OptionParser
from optparse import OptionGroup
import settings

def make_igvsession(igv_ed_umcu, igv_ed_hc, bam, vcf_hc, sample_id, vcf_SNV, axis, statistic):
    template_file = Template(open(settings.template_xml).read())
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
    parser = OptionParser()
    group = OptionGroup(parser, "Main options")
    group.add_option("--bam", dest = "bam", metavar = "[PATH]",
                     help = "bam file name"
                     )
    group.add_option("-o", dest = "output", metavar = "[PATH]",
                     help = "output_folder"
                     )
    group.add_option("--sampleid", dest = "sample_id", metavar = "[STRING]",
                     help = "sampleid name"
                     )
    group.add_option("--template", dest = "template", metavar = "[STRING]",
                     help = "Path to template XML"
                     )
    group.add_option("--refdate", dest = "ref_date", metavar = "[STRING]",
                     help = "Date of the used reference set)"
                     )
    group.add_option("--runid", dest = "run_id", metavar = "[STRING]",
                     help = "Name of run_id"
                     )
    parser.add_option_group(group)
    (opt, args) = parser.parse_args()

    output_folder = opt.output
    run_id = opt.run_id

    igv_ed_umcu = "UMCU/UMCU_{0}/UMCU_{1}_{0}_ref.igv".format(opt.bam, opt.ref_date)
    igv_ed_hc = "HC/HC_{0}/HC_{1}_{0}_ref.igv".format(opt.bam, opt.ref_date)
    bam_id = "../bam_files/{0}".format(opt.bam)
    vcf_hc = "HC/HC_{0}_{1}exome_calls.vcf".format(opt.ref_date, opt.bam)
    vcf_SNV = "../single_sample_vcf/{0}_{1}.vcf".format(opt.sample_id,run_id)
    igv_settings = settings.igv_settings
    for statistic in igv_settings:
        write_file = open("{0}/{1}_{2}_igv.xml".format(opt.output, opt.sample_id, statistic), "w")
        write_file.write(make_igvsession(igv_ed_umcu, igv_ed_hc, bam_id, vcf_hc, opt.sample_id, vcf_SNV, igv_settings[statistic], statistic))
        write_file.close()
