#! /usr/bin/env python3
import os
import sys
import subprocess
import re
import argparse
import glob
from multiprocessing import Pool
import pysam
import settings

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

def get_pu(bam, runid):
    """Get platform unit (PU) from bam file """
    merge = False
    workfile = pysam.AlignmentFile(bam, "rb")
    for readgroup in workfile.header['RG']:
        """ If one of the PU in all readgroups is not in the runID, the sample is considered to be a merge sample """
        if readgroup['PU'] not in runid:
             merge = True
    return merge

def multiprocess_ref(mp_list):


    action = "module load {renv} && Rscript {refscript} {folder}/ {folder}/{outputid} {targetbed} {refgenome} {exonbed}\n".format(
        renv = settings.r_version,
        refscript = settings.create_refset_r,
        folder = mp_list[0],
        outputid = mp_list[1],
        targetbed = mp_list[2],
        refgenome = settings.reference_genome,
        exonbed = mp_list[3]
        )
    os.system(action)

def make_refset(args):

    """Log all settings in setting.log file"""
    log_file="{output}/settings.log".format(
        output = args.output
        )
    write_file = open(log_file, "w")
    options = vars(args)
    for item in options:
        write_file.write("{0}\t{1}\n".format(str(item), str(options[item])))
    for item in dir(settings):
        if "__" not in item:
            write_file.write("{0}\t{1}\n".format(item, str(repr(eval("settings.%s" % item)))))
    write_file.close()


    """Make new reference set."""
    analysis = settings.analysis
    output_folder = format(os.path.abspath(args.output))

    bams = glob.glob("{}/**/*.bam".format(args.inputfolder), recursive=True)
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

    mp_list=[]
    for model in analysis:
        for item in ref_gender_dic:
            folder = "{0}/{1}_{2}_{3}".format(output_folder, model, item, str(args.output).rstrip("/").split("/")[-1])
            output_id = "{0}_{1}_{2}.EDref".format(model, item, args.prefix)
            action = "mkdir -p {0}".format(folder)
            os.system(action)
            for bam in ref_gender_dic[item]:
                os.system("ln -sd " + str(bam) + "* " + str(folder))
            mp_list+=[[folder, output_id, analysis[model]["target_bed"], analysis[model]["exon_bed"]]]

    with Pool(processes=int(args.simjobs)) as pool:
        result = pool.map(multiprocess_ref, mp_list, 1)

def multiprocess_call(multiprocess_list):

    """Log all settings in setting.log file"""
    log_file="{output}/{model}_{refset}_{bam}_{run}_settings.log".format(
        output = multiprocess_list[1],
        model = multiprocess_list[0],
        refset = args.refset,
        bam = args.sample,
        run = args.run
        )
    write_file = open(log_file, "w")
    options = vars(args)
    for item in options:
        write_file.write("{0}\t{1}\n".format(str(item), str(options[item])))
    for item in dir(settings):
        if "__" not in item:
            write_file.write("{0}\t{1}\n".format(item, str(repr(eval("settings.%s" % item)))))
    write_file.close()

    """Perform ExomeDepth analysis"""
    refset_R = "{refset_dir}/{model}_{gender}_{refset}.EDref".format(
        refset_dir = settings.refset_dir,
        model = multiprocess_list[0],
        gender = multiprocess_list[2],
        refset = args.refset
        )

    action = "module load {rversion} && Rscript {ed_r} {refset_R} {target_bed} {refgenome} {exon_bed} {prob} {bam} {model} {refset} {expected} {run}".format(
        rversion = settings.r_version,
        ed_r = settings.call_cnv_r,
        refset_R = refset_R,
        target_bed = multiprocess_list[3],
        refgenome = settings.reference_genome,
        exon_bed = multiprocess_list[4],
        prob = settings.probability[str(multiprocess_list[0])],
        bam =  multiprocess_list[5],
        model =  multiprocess_list[0],
        refset = args.refset,
        expected = args.expectedCNVlength,
        run = args.run 
        )
    os.system(action)

    """Perform csv to vcf conversion """
    action = "python {csv2vcf} {inputcsv} {refset} {model} {gender} {sampleid} {template}".format(
        csv2vcf = settings.csv2vcf,
        inputcsv = "{0}/{1}_{2}_{3}_{4}_exome_calls.csv".format(multiprocess_list[6],multiprocess_list[0],args.refset,args.inputbam,args.run),
        refset = args.refset,
        model = multiprocess_list[0],
        gender = multiprocess_list[2],
        sampleid = args.sample,
        template = settings.vcf_template
        )
    os.system(action)

def call_cnv(args):

    """Call CNV from BAMs"""
    bam = format(os.path.abspath(args.inputbam))
    output_folder = format(os.path.abspath(args.output))
    analysis = settings.analysis
    prob = settings.probability
    if not os.path.isdir(output_folder):
         os.system("mkdir -p " + str(output_folder)) 

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


    multiprocess_list=[]
    for model in analysis:
        multiprocess_list += [[model, output_folder, gender, analysis[model]["target_bed"], analysis[model]["exon_bed"], bam, output_folder]] 
    
    with Pool(processes=int(args.simjobs)) as pool:
        result = pool.map(multiprocess_call, multiprocess_list, 1)

    """Make log for stats of each model """
    merge = (get_pu(bam, args.run))
    for model in analysis:
        write_file = open("{output}/logs/{model}_{sample}_stats.log".format(output=args.output, model=model, sample=args.sample),"w")
        """ Get stats from VCF """
        vcf = "{output}/{model}/{model}_{refset}_{bam}_{run}_exome_calls.vcf".format(
            output=args.output,
            model=model,
            refset=args.refset,
            bam=bam.split("/")[-1],
            run=args.run
            )
        stats = (subprocess.getoutput("tail -n1 {}".format(vcf)).split()[-1]).split(":")
        correlation, del_dup_ratio, number_calls = float(stats[4]), float(stats[8]), float(stats[9])

        printline = ("{sample}\t{model}\t{correlation}\t{del_dup_ratio}\t{number_calls}".format(
            sample=args.sample,
            model=model,
            correlation=correlation,
            del_dup_ratio=del_dup_ratio,
            number_calls=int(number_calls)
            ))

        if args.qc_stats:
            failed = False
            if correlation < settings.correlation or number_calls > settings.number_calls or del_dup_ratio < settings.del_dup_ratio[0] or del_dup_ratio > settings.del_dup_ratio[1]:
                failed = True
                printline += "\tFAIL"
            else:
                printline += "\tOK"

        write_file.write(printline + "\n")
        write_file.close()    


    """Make IGV session xml """
    action = "python {igv_xml} {bam} {output} {sampleid} {template} {refdate} {runid} --pipeline {pipeline}".format(
        igv_xml = settings.igv_xml,
        bam = args.inputbam,
        output = args.output,
        sampleid = args.sample,
        template = settings.template_xml,
        refdate = args.refset,
        runid = args.run,
        pipeline = args.pipeline
        )
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
    parser_refset.add_argument('--simjobs', default=4, help='number of simultaneous samples to proces. Note: make sure similar threads are reseved in session! [default = 4]')
    parser_refset.add_argument('--genderfile', help='Gender file: tab delimited txt file with bam_id  and gender (as male/female)')
    parser_refset.set_defaults(func = make_refset)

    parser_cnv = subparser.add_parser('callcnv', help='Call CNV with ExomeDepth basedon BAM file')
    parser_cnv.add_argument('output', help='output folder for CNV calling')
    parser_cnv.add_argument('inputbam', help='Input BAM file')
    parser_cnv.add_argument('run', help='Name of the run')
    parser_cnv.add_argument('sample', help='Sample name')
    parser_cnv.add_argument('refset', help='Reference set to be used (e.g. Jan2020)')
    parser_cnv.add_argument('--pipeline', default='nf', choices=['nf', 'iap'], help='pipeline used for sample processing (nf = nexflow (default), IAP = illumina analysis pipeline')
    parser_cnv.add_argument('--simjobs', default=2, help='number of simultaneous samples to proces. Note: make sure similar threads are reseved in session! [default = 2]')
    parser_cnv.add_argument('--genderfile', help='Gender file: tab delimited txt file with bam_id  and gender (as male/female)')
    parser_cnv.add_argument('--expectedCNVlength',default=settings.expectedCNVlength, help='expected CNV length (basepairs) taken into account by ExomeDepth [default expectedCNVlength in settings.py]')
    parser_cnv.add_argument('--qc_stats', action='store_true', help='switch on QC check for exomedepth VCF stats (default = off)')
    parser_cnv.set_defaults(func = call_cnv)

    args = parser.parse_args()
    args.func(args)

