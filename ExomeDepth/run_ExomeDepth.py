#! /usr/bin/env python3
import os
import sys
import pysam
import subprocess
import re
import settings
import argparse

def valid_read(read):
    """Check if a read is properly mapped."""
    if (read.mapping_quality >= 20 and read.reference_end and read.reference_start):
        return True
    else:
        return False


def get_gender(bam):
    """Determine chrY ratio based on read count in bam (excl PAR)."""
    workfile = pysam.AlignmentFile(bam, "rb")
    locus = 'Y:2649520-59034050'
    yreads = float(sum([valid_read(read) for read in workfile.fetch(region=locus)]))
    total = float(workfile.mapped)
    yratio = float("%.2f" % ((yreads / total) * 100))
    if yratio <= 0.06:
        return "female"
    elif yratio >= 0.12:
        return "male"
    else:
        return "unknown"

def make_refset(args):

    """Make new reference set."""
    analysis = settings.analysis
    bams = subprocess.getoutput("find -L {0} -iname \"*.realigned.bam\"".format(args.inputfolder)).split()
    print("Number of BAM files detected = {0}".format(len(bams)))
 
    """Get gender from chrY read count ratio."""
    ref_gender_dic = {}  #Dictionary with gender of each sample
    for bam in bams:
       gender = get_gender(bam)
       if gender is not "unknown":
           if gender not in ref_gender_dic:
               ref_gender_dic[gender] = [bam]      
           else:
               ref_gender_dic[gender] += [bam]
       else:
           print("Sample {0} has unknown gender and is removed from analysis".format(bam))
   
    """Make folder per gender + analysis, and soflink BAMs in these folders."""
    for model in analysis:
        for item in ref_gender_dic:
            folder = "{0}/{1}_{2}_{3}".format(args.output, model, item, str(args.output).split("/")[-1])
            output_id = "{0}_{1}_{2}.EDref".format(model, item, args.prefix)
            action = "mkdir -p {0}".format(folder)
            os.system(action)
            for bam in ref_gender_dic[item]:
                os.system("ln -sd " + str(bam) + "* " + str(folder))
            action = ("module load {0} && Rscript {1} {2}/ {2}/{3} {4} {5} {6}\n".format(
                settings.r_version,
                settings.create_refset_r,
                folder,
                output_id,
                analysis[model]["target_bed"],
                settings.reference_genome,
                analysis[model]["exon_bed"]
                ))
            os.sytem(action)

def call_cnv(args):

    """Call CNV from BAMs"""
    bam = format(os.path.abspath(args.inputbam))
    output_folder = format(os.path.abspath(args.output))
    analysis = settings.analysis
    ed_r = settings.call_cnv_r
    refgenome = settings.reference_genome
    prob = settings.probability
    
    """Determine gender"""
    gender = get_gender(bam)
    if args.genderfile:  # overrule gender as given in gender_file
       gender_dic = gender_file(args.genderfile)
       for item in gender_dic:
            if str(item) in str(bam):
               gender = gender_dic[item]
    if gender == "unknown":  # try to determine gender on BAM ID annotation
        if re.search('[C|P]M', bam.split("/")[-1]):
            print("Sample {0} has a unknown gender based on chrY reads, but resolved as male based on sampleID".format(bam.split("/")[-1]))
            gender = "male"
        elif re.search('[C|P]F', bam.split("/")[-1]):
            print("Sample {0} has a unknown gender based on chrY reads, but resolved as female based on sampleID".format(bam.split("/")[-1]))
            gender = "female"
        else:
            sys.exit("Sample {0} has a unknown gender and will not be analysed".format(bam.split("/")[-1]))

    for model in analysis:
        """Log all settings in setting.log file"""
        log_file="{0}_{1}_{2}_settings.log".format(model,args.refset,args.inputbam)
        write_file = open(log_file, "w")
        options = vars(args) 
        for item in options:
            write_file.write("{0}\t{1}\n".format(str(item), str(options[item])))
        for item in dir(settings):
            if "__" not in item:
                write_file.write("{0}\t{1}\n".format(item, str(repr(eval("settings.%s" % item)))))
        write_file.close()

        """Perform ExomeDepth analysis"""
        action = ("module load {0} && Rscript {1} {2} {3} {4} {5} {6} {7} {8} {9}".format(
            settings.r_version,
            ed_r,
            analysis[model]["refset"][gender],
            analysis[model]["target_bed"],
            refgenome,
            analysis[model]["exon_bed"],
            prob[str(model)],
            bam,
            model,
            args.refset
            ))
        os.system(action)

        """Perform csv to vcf conversion """
        #csv_file="{0}/{1}_{2}_{3}_exome_calls.csv".format(output_folder,model,args.refset,args.inputbam)
        action = ("python {csv2vcf} {inputcsv} {refset} {model} {gender} {sampleid} {template}".format(
            csv2vcf = settings.csv2vcf,
            inputcsv = "{0}/{1}_{2}_{3}_exome_calls.csv".format(output_folder,model,args.refset,args.inputbam),
            refset = args.refset,
            model = model,
            gender = gender,
            sampleid = args.sample,
            template = settings.vcf_template
            ))
        os.system(action)

    """Make IGV session xml """
    action = ("python {igv_xml} {bam} {output} {sampleid} {template} {refdate} {runid}".format(
        igv_xml = settings.igv_xml,
        bam = args.inputbam,
        output = output_folder,
        sampleid = args.sample,
        template = settings.igv_xml,
        refdate = args.refset,
        runid = args.run
        ))
    os.system(action)

def gender_file(genderfile):
    gender_dic = {}
    gender_file = open(str(genderfile), "r")
    for line in gender_file:
        splitline = line.split()
        if splitline[0] not in gender_dic:
            gender_dic[splitline[0]] = splitline[1]
    gender_file.close()
    return gender_dic

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()
    
    parser_refset = subparser.add_parser('makeref', help='Make ExomeDepth reference set')
    parser_refset.add_argument('output', help='Output folder for reference set files')
    parser_refset.add_argument('inputfolder', help='Input folder containing BAM files')
    parser_refset.add_argument('prefix', help='Prefix for reference set (e.g. Jan2020)')
    parser_refset.add_argument('--genderfile', help='Gender file: tab delimited txt file with bam_id  and gender (as male/female)')
    parser_refset.set_defaults(func = make_refset)

    parser_cnv = subparser.add_parser('callcnv', help='Call CNV with ExomeDepth basedon BAM file')
    parser_cnv.add_argument('output', help='output folder for CNV calling')
    parser_cnv.add_argument('inputbam', help='Input BAM file')
    parser_cnv.add_argument('run', help='Name of the run')
    parser_cnv.add_argument('sample', help='Sample name')
    parser_cnv.add_argument('refset', help='Reference set to be used (e.g. Jan2020)')
    parser_cnv.add_argument('--genderfile', help='Gender file: tab delimited txt file with bam_id  and gender (as male/female)')
    parser_cnv.set_defaults(func = call_cnv)

    args = parser.parse_args()
    args.func(args)

