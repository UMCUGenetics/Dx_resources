#! /usr/bin/env python
import sys
import os
from string import Template
from optparse import OptionParser
from optparse import OptionGroup
import settings

def make_igvsession(igv_ed_umcu, igv_ed_hc, bam, vcf_hc, sample_id, vcf_SNV):
    template_file = Template(open(settings.template_xml).read())
    new_session = str(sample_id) + "_igv.xml"
    igv_ed_hc_ratio_test = str(igv_ed_hc) + "_ratio_test"
    igv_ed_umcu_ratio_test = str(igv_ed_umcu) + "_ratio_test"
    bam_coverage = str(bam) + "_coverage"
    bam_junctions = str(bam) + "_junctions"
    substitute_dic={'session_var' : new_session, 'igv_ed_umcu' : igv_ed_umcu, 
                    'igv_ed_hc' : igv_ed_hc, 'bam':bam, 'vcf_hc' : vcf_hc, 
                    'sample_id' : sample_id, 'igv_ed_hc_ratio_test' : igv_ed_hc_ratio_test,
                    'igv_ed_umcu_ratio_test' : igv_ed_umcu_ratio_test, 
                    'bam_coverage' : bam_coverage, 'bam_junctions' : bam_junctions, 
                    'vcf_SNV': vcf_SNV
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
    write_file = open("{0}/{1}_igv.xml".format(opt.output,opt.sample_id), "w")
    write_file.write(make_igvsession(igv_ed_umcu, igv_ed_hc, bam_id, vcf_hc, opt.sample_id, vcf_SNV))
    write_file.close()
