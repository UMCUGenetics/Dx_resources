#! /usr/bin/env python
import sys
import os
import re
import subprocess
from optparse import OptionParser
from optparse import OptionGroup
import vcf
import collections
import copy
import pysam
import pandas as pd
import settings
import decimal

def cnv_locationtype(region, par1, par2):
    chrom = str(region[0]).upper()
    start = int(region[1])
    stop = int(region[2])
    if chrom == "X":
        if start >= par1[0] and stop <= par1[1] or start >= par2[0] and stop <= par2[1]:  #If CNV is nested in par1 or par2 region
            return "chrXpar"
        else:
            return "chrX"
    elif chrom == "Y":
        return "chrY"
    else:
        return "auto"


if __name__ == "__main__":
    parser = OptionParser()
    group = OptionGroup(parser, "Main options")
    group.add_option("--csv", dest = "input_csv", metavar = "[PATH]",
                     help = "full path to input CSV file"
                     )
 
    group.add_option("--template", dest = "template", metavar = "[STRING]",
                     help = "Path to template VCF"
                     )

    group.add_option("--refset", dest = "refset", metavar = "[STRING]",
                     help = "used refset"
                     )

    group.add_option("--model", dest = "model", metavar = "[STRING]",
                     help = "used model"
                     )
    group.add_option("--gender", dest = "gender", metavar = "[STRING]",
                     help = "used gender"
                     )
    group.add_option("--sampleid", dest = "sampleid", metavar = "[STRING]",
                     help = "sampleid name to be included in VCF"
                     )
    parser.add_option_group(group)
    (opt, args) = parser.parse_args()


    if not opt.input_csv:
         sys.exit("provide input file (--csv)")

    if opt.template:
        if os.path.isfile(str(opt.template)):
            template_vcf = str(opt.template)
        else:
            sys.exit("Input template VCF does not exist")
    else:
        sys.exit("provide template VCF (-t)")

    sampleid= opt.sampleid

    vcf_reader = vcf.Reader(open(template_vcf, 'r'))
    format_keys = vcf_reader.formats.keys()
    record = vcf_reader.__next__()  # First record (dummy for all other records)
    new_record = copy.deepcopy(record)
    new_record.samples[0].data = collections.namedtuple('CallData', format_keys)  # For single sample VCF only!
    format_vals = [record.samples[0].data[vx] for vx in range(len(format_keys))]
    format_dict = dict(zip(format_keys, format_vals))
    for f in ['GT', 'CCN', 'BF', 'RT', 'CR', 'RS', 'IH', 'CM', 'PD', 'TC']:
        format_dict[f] = ""
    new_vals = [format_dict[x] for x in format_keys]
    new_record.samples[0].data = new_record.samples[0].data._make(new_vals)

    df_csv = pd.read_csv(opt.input_csv)
    vcf_reader.samples = [sampleid]  # Change template sampleid in sampleid
    """Add metadata ED reference set used."""
    vcf_reader.metadata['EDreference'] = ["{0}_{1}_{2}".format(opt.model,opt.gender,opt.refset)]
    with open(opt.input_csv[0:-4]+".vcf", 'w') as vcf_output_file:
        vcf_writer = vcf.Writer(vcf_output_file, vcf_reader)

        """Determine percentage DEL/(DEL+DUP) for all calls in VCF."""
        dels = 0
        dups = 0
        for index, row in df_csv.iterrows():
            if row['type'] == "deletion":
                dels += 1
            if row['type'] == "duplication":
                 dups += 1
        perc_del = "%.2f" % ((float(dels) / (float(dels) + float(dups))) * 100)
        total_calls = dels + dups

        for index, row in df_csv.iterrows():  # index not used as we only use single sample VCF
            """Change record fields."""
            new_record.CHROM = row['chromosome']
            new_record.POS = row['start']
            new_record.ID = "."
            row_type = str(row['type'])

            """Include reference genome base"""
            ref_file = vcf_reader.metadata['reference']
            fasta_open = pysam.Fastafile(ref_file)
            seq_fasta = fasta_open.fetch(str(row['chromosome']), int(row['start']-1), int(row['start']))  # 0-based coordinates
            new_record.REF = seq_fasta
 
            """Write type of call."""
            if row_type == "duplication":
                new_record.ALT = ["<DUP>"]
                new_record.INFO['SVTYPE'] = "DUP"
            elif row_type == "deletion":
                new_record.ALT = ["<DEL>"]
                new_record.INFO['SVTYPE'] = "DEL"
            else:
                new_record.ALT = ["NaN"]
                new_record.INFO['SVTYPE'] = "NaN"
            new_record.QUAL = "1000"  # as STRING
            new_record.FILTER = "PASS"

            """Determine genotype."""
            ratio = row['reads.ratio']
            if str(ratio) == "inf":  #Rename infinitity values to 99
                ratio = 99

            """Consider homozygous genotype only for deletion and with ratio <0.25."""
            if str(row['type']).lower() == "deletion" and float(ratio) < float(settings.ratio_threshold_del):
                genotype = "1/1"
            else:  # Always het for duplication, and het for deletion if not < settings.ratio_threshold_del
                genotype = "0/1"

            """Determine copy number. Note this will not work for mosaik events"""
            par1 = settings.par1
            par2 = settings.par2
            normal_CN = settings.normal_CN
            region = [str(row['chromosome']), int(row['start']), int(row['end'])]
            locus_type = cnv_locationtype(region, par1, par2)
            normal_copy = float(normal_CN[opt.gender][locus_type])
            calc_copynumber = normal_copy  * float(ratio)

            # Estimate true copynumber by rounding to nearest integer
            copynumber = int(decimal.Decimal(calc_copynumber).quantize(decimal.Decimal('0'), rounding=decimal.ROUND_HALF_UP))
            if opt.gender == "female" and locus_type == "chrY":
                """In case CNV is detected on chrY in female, correct for this"""
                print ("WARNING: {sample} chromosome Y CNV detected (region = {region}) in female, calc_copynumber set to 0 (deletion call) or 1 (duplication call)".format(
                       sample = str(sampleid),
                       region = str("_".join(str(x) for x in region))
                      ))    

                # CNV CN is set to 1, could also be >1 (but makes no biological sense)
                if ratio > 1:
                    calc_copynumber = 1
                    genotype = "1/1"
                else:  # Assuming CNV is called by noise/mismapping/ on chrY, set CN to 0. Could result in masking contamination of male sample?
                    calc_copynumber = 0
                    genotype = "1/1"
            else:
                if row_type == "deletion" and calc_copynumber > normal_copy or row_type == "duplication" and calc_copynumber < normal_copy:
                    """ If calc_copynumber is opposite of expected CN for region, i.e. ratio 1.5 for a deletion"""
                    print ("WARNING: {sample} CNV copynumber estimation {copynumber} does not match CNV type {rowtype} for region {region}".format(
                           sample = str(sampleid),
                           copynumber = str(float(calc_copynumber)),
                           rowtype = row_type,
                           region =  str("_".join(str(x) for x in region))
                           ))
                    """Note: no correction here. should be bugfix in the ExomeDepth code"""
  
                if copynumber == int(normal_copy):
                    """ Estimated copynumber is similar to copyneutral """
                    print ("WARNING: {sample} true copynumber for region {region} is same as normal CN > set to -1 for deletion, +1 for duplication".format(
                           sample = str(sampleid),
                           region =  str("_".join(str(x) for x in region))
                          ))
                    if row_type == "deletion":  # If deletion correct copynumber with -1
                        copynumber -= 1
                        if copynumber == 0:
                            genotype = "1/1"
                    elif row_type == "duplication":  #If duplication correct copynumber with +1
                        copynumber += 1
  
            """Change INFO fields"""
            new_record.INFO['END'] = row['end']
            new_record.INFO['NTARGETS'] = row['nexons']
            new_record.INFO['SVLEN'] = int(row['end']) - int(row['start'])  # Input is assumed 0-based
            new_record.INFO['CN'] = copynumber
            call_conrad = row['Conrad.hg19']
            if str(call_conrad) == "nan":
                call_conrad = "NaN"
            new_record.INFO['cCNV'] = call_conrad

            """Change FORMAT fields"""
            for f in ['GT', 'CCN', 'BF', 'RT', 'CR', 'RS', 'IH', 'CM', 'PD', 'TC']:
                format_dict[f] = ""
            format_dict['GT'] = str(genotype)
            format_dict['CCN'] = "%.2f" % (float(calc_copynumber))
            format_dict['BF'] = "%.2f" % (float(row['BF']))
            format_dict['RT'] = "%.2f" % (float(ratio))
            format_dict['CR'] = "%.4f" % (float(row['correlation']))
            format_dict['RS'] = row['refsize']
            format_dict['IH'] = "NaN"  # Inheritence is not build in yet
            format_dict['CM'] = opt.model
            format_dict['PD'] = perc_del
            format_dict['TC'] = total_calls
            new_vals = [format_dict[x] for x in format_keys]
            new_record.samples[0].data = new_record.samples[0].data._make(new_vals)  # NOTE: GT must be first in order of metadata!
            vcf_writer.write_record(new_record)

        vcf_writer.flush()
