#! /usr/bin/env python
import sys, os,re, commands
from optparse import OptionParser
from optparse import OptionGroup
import pathlib
import vcf
import collections
import copy
import pysam
import numpy as np
import pandas as pd
import settings

def CNV_LocationType(region):
    if region[0]=="X":
        if region[1] >= PAR1[0] and region[2] <= PAR1[1] or region[1] >= PAR2[0] and region[2] <= PAR2[1]:      # if CNV is nested in PAR1 or PAR2 region!
            return "chrXPAR"
        else:
            return "chrX"

    elif region[0]=="Y":
        return "chrY"
    else:
        return "auto"

if __name__ == "__main__":
        parser = OptionParser();
        group = OptionGroup(parser, "Main options")
        group.add_option("-i", dest="input", metavar="[PATH]", help="input folder including csv files")
        group.add_option("-t",dest="template", metavar="[STRING]", help="Path to template VCF")
        group.add_option("-m",dest="callmodel", metavar="[STRING]", help="Name of calling model used (order: [Model]_[Gender]_[Date].EDRef)")
        parser.add_option_group(group)
        (opt, args) = parser.parse_args()

if opt.input:
    csv_dir =str(opt.input)
else:
    sys.exit("provide input folder (-i)")

if opt.template:
    if os.path.isfile(str(opt.template)):
        template_vcf=str(opt.template)
    else:
        sys.exit("Input template VCF does not exist")
else:
    sys.exit("provide template VCF (-t)")

# determine model, gender, date
if opt.callmodel:
    referenceset=str(opt.callmodel).split("/")[-1]
    model=referenceset.split("_")[0]
    gender=referenceset.split("_")[1]
else:
    #model=str(csv.split("/")[-2].split("_")[0])                                     # Model name based on output folder (might not be stable)
    sys.exit("please provide referenceset used \"-m\" ")


csvs=commands.getoutput("find -L "+str(csv_dir)+" -iname \"*.csv\"").split()

vcf_reader = vcf.Reader(open(template_vcf, 'r'))                                    	# load template VCF
f_keys = vcf_reader.formats.keys()                                                  	# it's an ordered dict
record = vcf_reader.next()                                                          	# First record (dummy for all other records)
new_record = copy.deepcopy(record)                                                  	# make new record
new_record.samples[0].data = collections.namedtuple('CallData', f_keys)			# For single sample VCF only!
f_vals = [record.samples[0].data[vx] for vx in range(len(f_keys))]
handy_dict = dict(zip(f_keys, f_vals))
for f in ['GT','CCN','BF','RT','CR','RS','IH','CM','PD','TC'] :
    handy_dict[f] = ""                                                               # String for modification ("empty")
new_vals = [handy_dict[x] for x in f_keys]
new_record.samples[0].data = new_record.samples[0].data._make(new_vals)

for csv in csvs:
    csv_file=open(csv,"r")
    sampleID=csv.split("/")[-1].split("_")[0] 
    df_csv = pd.read_csv(csv)
    vcf_reader.samples=[sampleID]							# set correct sample name
    vcf_reader.metadata['EDreference']= [str(opt.callmodel)]				# add metadata ED reference set used
    vcf_writer = vcf.Writer(open(str(csv)[0:-4]+".vcf", 'w'),vcf_reader)             	# write output VCF

    # determine percentage DEL/(DEL+DUP) for all calls in VCF
    dels=0
    dups=0
    for index, row in df_csv.iterrows():
        if row['type']== "deletion":
            dels+=1 
        if row['type']== "duplication":
            dups+=1
    perc_del="%.2f" % ((float(dels)/(float(dels)+float(dups)))*100)
    total_calls= dels+dups

    for index, row in df_csv.iterrows():						# index not used as we only use single sample VCF
        ## Change record fields
        new_record.CHROM=row['chromosome']
        new_record.POS=row['start']
        new_record.ID="."
        
        # include reference genome base here
        ref_file= vcf_reader.metadata['reference']
        fasta_open = pysam.Fastafile(ref_file)
        seq_fasta = fasta_open.fetch(str(row['chromosome']), int(row['start']-1), int(row['start']))	# 0-based coordinates
        new_record.REF=seq_fasta						

        if row['type']== "duplication":							# <DUP> or <DEL>
            new_record.ALT=["<DUP>"]							
            new_record.INFO['SVTYPE']="DUP"
        elif row['type']== "deletion":
            new_record.ALT=["<DEL>"]
            new_record.INFO['SVTYPE']="DEL"
        else:
            new_record.ALT=["NaN"]
            new_record.INFO['SVTYPE']="NaN"
        new_record.QUAL="1000"
        new_record.FILTER="PASS"
 
        # determine genotype
        ratio=row['reads.ratio']
        if str(ratio) == "inf":                                                         # rename infinitity values to 99
            ratio=99
        if float(ratio) < float(settings.ratio_threshold_del) and str(row['type'])== "deletion":       # Consider homozygous genotype only for deletion and ratio <0.25
            genotype="1/1"
        else:                                                                           # Otherwise always heterozygous.
            genotype="0/1"

        # determine copy number. ### NOT this will not work for mosaik events
        PAR1=settings.PAR1
        PAR2=settings.PAR2
        normal_CN=settings.normal_CN
        region=[str(row['chromosome']),int(row['start']),int(row['end'])]
        locus_type= CNV_LocationType(region)
        calc_copynumber=float(normal_CN[gender][locus_type])*float(ratio) 
        copynumber=int(round(calc_copynumber)) 						# Estimate true copynumber by rounding to nearest integer
        if gender == "female" and locus_type == "chrY":
            print "WARNING: "+str(sampleID)+" chromosomeY CNV detected (region = "+str("_".join(str(x) for x in region))+") in female, calc_copynumber set to 0 (deletion call) or 1 (duplication call)"		# this is not correct. need to resolve this (?)
            if ratio > 1:
                calc_copynumber=1
                genotype="1/1"
            else:
                calc_copynumber=0
                genotype="1/1"
        else:
            if str(row['type'])== "deletion" and calc_copynumber > float(normal_CN[gender][locus_type]) or str(row['type'])== "duplication" and calc_copynumber < float(normal_CN[gender][locus_type]):	#no bugfixed included yet
                print "WARNING: "+str(sampleID)+" CNV copynumber estimation "+str(float(calc_copynumber))+" does not match CNV type "+str(row['type'] + " for region "+str("_".join(str(x) for x in region)))

            if copynumber == int(normal_CN[gender][locus_type]):
                print "WARNING, "+str(sampleID)+" true copynumber for region "+str("_".join(str(x) for x in region))+" is same as normal CN > set to -1 for deletion, +1 for duplication"
                if str(row['type'])== "deletion":
                    copynumber-=1
                    if copynumber == 0:
                        genotype="1/1" 
                elif str(row['type'])== "duplication":
                    copynumber+=1

        ## Change INFO fields 
        new_record.INFO['END']=row['end'] 
        new_record.INFO['NEXONS']=row['nexons']
        new_record.INFO['SVLEN']=int(row['end'])-int(row['start'])                      # input is assumed 0-based
        new_record.INFO['CN']= copynumber
        call_conrad=row['Conrad.hg19']
        if str(call_conrad) == "nan":							# rename missing value to NaN
            call_conrad = "NaN"
        new_record.INFO['cCNV']=call_conrad
     
        ## Change FORMAT field
        for f in ['GT','CCN','BF','RT','CR','RS','IH','CM','PD','TC'] :
            handy_dict[f] = ""
        handy_dict['GT'] = str(genotype)
        #handy_dict['CN'] = copynumber
        handy_dict['CCN'] = calc_copynumber
        handy_dict['BF'] = "%.2f" % (float(row['BF']))
        handy_dict['RT'] = "%.2f" % (float(ratio))
        handy_dict['CR'] = "%.4f" % (float(row['correlation']))
        handy_dict['RS'] = row['refsize']
        handy_dict['IH'] = "NaN"								# inheritence is not build in yet
        handy_dict['CM'] = model
        handy_dict['PD'] = perc_del
        handy_dict['TC'] = total_calls
        new_vals = [handy_dict[x] for x in f_keys]
        new_record.samples[0].data = new_record.samples[0].data._make(new_vals)	# NOTE: GT must be first in order of metadata!
        vcf_writer.write_record(new_record)
        
    vcf_writer.flush()
