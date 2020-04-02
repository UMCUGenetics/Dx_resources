#! /usr/bin/env python
import sys
import os
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
    min_axis = axis[0]
    mid_axis = axis[1]
    max_axis = axis[2]

    substitute_dic={'session_var' : new_session, 'igv_ed_umcu' : igv_ed_umcu, 
                    'igv_ed_hc' : igv_ed_hc, 'bam':bam, 'vcf_hc' : vcf_hc, 
                    'sample_id' : sample_id, 'igv_ed_hc_test' : igv_ed_hc_test,
                    'igv_ed_umcu_test' : igv_ed_umcu_test, 
                    'bam_coverage' : bam_coverage, 'bam_junctions' : bam_junctions, 
                    'vcf_SNV': vcf_SNV, 'min_axis' : min_axis, 'mid_axis' : mid_axis,'max_axis' : max_axis
                   }
    new_file = template_file.substitute(substitute_dic)
    return new_file

if __name__ == "__main__":
    parser = OptionParser()
    group = OptionGroup(parser, "Main options")
    group.add_option("-b", dest = "bam", metavar = "[PATH]",
                     help = "bam file name"
                     )
    group.add_option("-o", dest = "output", metavar = "[PATH]",
                     help = "output_folder"
                     )
    group.add_option("-i", dest = "sample_id", metavar = "[STRING]",
                     help = "sampleid name"
                     )
    group.add_option("-t", dest = "template", metavar = "[STRING]",
                     help = "Path to template XML"
                     )
    group.add_option("-m", dest = "callmodel", metavar = "[STRING]",
                     help = "Name of calling model used (order: [Model]_[Gender]_[Date].EDRef)"
                     )
    parser.add_option_group(group)
    (opt, args) = parser.parse_args()

    output_folder = opt.output
    refdate = opt.callmodel.split("/")[-1].split("_")[2].split(".")[0]
    igv_ed_umcu = "UMCU/UMCU_{0}/UMCU_{1}_{0}_ref.igv".format(opt.bam, refdate)
    igv_ed_hc = "HC/HC_{0}/HC_{1}_{0}_ref.igv".format(opt.bam, refdate)
    bam_id = "../{0}/mapping/{1}".format(opt.sample_id,opt.bam)
    vcf_hc = "HC/HC_{0}_{1}exome_calls.vcf".format(refdate, opt.bam)
    vcf_SNV = "../fingerprint/{0}_fingerprint.vcf".format(opt.sample_id)
    
    igv_settings = settings.igv_settings
    for statistic in igv_settings:
        write_file = open("{0}/{1}_{2}_igv.xml".format(opt.output, opt.sample_id, statistic), "w")
        write_file.write(make_igvsession(igv_ed_umcu, igv_ed_hc, bam_id, vcf_hc, opt.sample_id, vcf_SNV, igv_settings[statistic], statistic))
        write_file.close()
