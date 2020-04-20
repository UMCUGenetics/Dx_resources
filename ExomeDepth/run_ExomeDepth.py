#! /usr/bin/env python
import os
import sys
import pysam
import subprocess
import re
from optparse import OptionParser, OptionGroup
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


if __name__ == "__main__":
    parser = OptionParser()
    group = OptionGroup(parser, "Main options")
    group.add_option("-r", dest="make_ref",
                     action = 'store_true', help="Make new reference set")
    group.add_option("-c", dest = "make_call",
                     action = 'store_true', help = "Call CNV from BAMs")

    group.add_option("-o", dest = "output_folder",
                     default = "./", metavar = "[PATH]",
                     help = "path to output folder [default = ./]")
    group.add_option("--bam", dest = "input_bam", 
                     metavar = "[PATH]", help = "path to input bam file")
    group.add_option("--run", dest = "run_id",
                     metavar = "[STRING]", help = "run_id name")
    group.add_option("--model", dest = "model", default = "HC",
                     metavar = "[STRING]", help = "model HC/UMCU [default = HC]")
    group.add_option("--sample", dest = "sample_id", 
                     metavar = "[STRING]", help = "sampleID")


    group.add_option("--if", dest = "input_folder",
                     metavar = "[PATH]", help = "folder path in which contains BAM files")
    group.add_option("-g", dest = "gender_file",
                     metavar = "[PATH]",
                     help = "full path to tab delimited txt file with bam \
                    + gender (male/female) [default = off]")
    group.add_option("-p", dest = "prefix", metavar = "[PATH]",
                     help = "referenceset naming [default = input_folder basename]")
    group.add_option("--refset", default= settings.refset, dest = "refset", metavar = "[STRING]",
                     help = "reference set to be used [default = reference set in setting.py]")
    parser.add_option_group(group)
    (opt, args) = parser.parse_args()
   

    """ 
    if opt.input_bam and not opt.make_ref:
        wkdir = "/".join(opt.input_bam.split("/")[0-1]).rstrip("/")
    elif opt.input_folder:
        wkdir = opt.input_folder.rstrip("/") 
    else:
        sys.exit("please provide input folder/bam")

    if not opt.output_folder:
        outdir = format(os.path.abspath("./"))    
    else:
        outdir = format(os.path.abspath(opt.output_folder))

    if not os.path.isdir(outdir):
        os.system("mkdir -p " + str(outdir))
    
    if opt.prefix:
        prefix = opt.prefix
    else:
        prefix = str(wkdir).split("/")[-1]

    if opt.run_id:
        run_id=opt.run_id  
   
    if opt.sample_id:
        sample_id= opt.sample_id
    else:
        sys.exit("please provide sample_id")
    """

    bam = format(os.path.abspath(opt.input_bam))
 
    if not opt.output_folder:
        output_folder = format(os.path.abspath("./"))
    else:
        output_folder = format(os.path.abspath(opt.output_folder))

    output_folder_model = "{0}/{1}".format(output_folder,opt.model)
    if not os.path.isdir(output_folder_model):
        os.system("mkdir "+ output_folder_model)

    analysis = settings.analysis
    ed_r = settings.call_cnv_r
    refgenome = settings.reference_genome
    prob = settings.probability
    gender = settings.gender

    if opt.gender_file:
        gender_dic = {}
        gender_file = open(str(opt.gender_file), "r")
        for line in gender_file:
            splitline = line.split()
            if splitline[0] not in gender_dic:
                gender_dic[splitline[0]] = splitline[1]
        gender_file.close()
   
    """Log all settings in setting.log file"""
    if opt.input_bam:
        log_file="{0}_{1}_{2}_settings.log".format(opt.model,opt.refset,opt.input_bam)
        write_file = open(log_file, "w") 
        (options, args) = parser.parse_args()
        for item in vars(options):
            write_file.write("{0}\t{1}\n".format(str(item), str(vars(options)[item])))
        for item in dir(settings):
            if "__" not in item:
                write_file.write("{0}\t{1}\n".format(item, str(repr(eval("settings.%s" % item)))))
        write_file.close()
  
    if opt.make_ref and not opt.make_call:  
        """Make new reference set."""
        bams = subprocess.getoutput("find -L {0} -iname \"*.realigned.bam\"".format(wkdir)).split()
        print("Number of BAM files detected = {0}".format(len(bams)))
        """Get gender from chrY read count ratio."""

        ref_gender_dic = {}  #Dictionary with gender of each sample
        for item in gender:
            if str(item) not in ref_gender_dic:
                ref_gender_dic[str(item)] = []
    
        for bam in bams:
            gender = get_gender(bam)
            if gender is not "unknown":
                ref_gender_dic[get_gender(bam)] += [bam]
            else:
                print("Sample {0} has unknown gender and is removed from analysis".format(bam))
        """Make folder per gender + analysis, and soflink BAMs in these folders."""
        for target in analysis:
            for item in ref_gender_dic:
                folder = "{0}/{1}_{2}_{3}".format(output_folder, target, item, str(wkdir).split("/")[-1])  
                output_id = "{0}_{1}_{2}.EDref".format(target, item, prefix)
                os.system("mkdir -p {0}".format(folder))
                for bam in ref_gender_dic[item]:
                    os.system("ln -sd " + str(bam) + "* " + str(folder))
                os.system("module load {0} && Rscript {1} {2}/ {2}/{3} {4} {5} {6}\n".format(
                                 settings.r_version,
                                 settings.create_refset_r,
                                 folder,
                                 output_id,
                                 analysis[target]["target_bed"],
                                 settings.reference_genome,
                                 analysis[target]["exon_bed"]
                                 ))

    elif opt.make_call and not opt.make_ref:  
        """Call CNV from BAMs"""
        gender = get_gender(bam)
        if opt.gender_file:  # overrule gender as given in gender_file
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

        """Perform ExomeDepth analysis"""
        action = ("module load {0} && Rscript {1} {2} {3} {4} {5} {6} {7} {8} {9}".format(
                   settings.r_version,
                   ed_r,
                   analysis[opt.model]["refset"][gender],
                   analysis[opt.model]["target_bed"],
                   refgenome,
                   analysis[opt.model]["exon_bed"],
                   prob[str(opt.model)],
                   bam,
                   opt.model,
                   opt.refset
                  ))

        print(action)
        os.system(action)

        """Move output to specifiek folder"""
        os.system("mv {0}/* {1}".format(output_folder,output_folder_model))
        
        """Perform csv to vcf conversion """
        csv_file="{0}/{1}_{2}_{3}_exome_calls.csv".format(output_folder_model,opt.model,opt.refset,opt.input_bam)
        action = ("python {0} --csv {1} --refset {2} --model {3} --gender {4} --sampleid {5} --template {6}".format(
                   settings.csv2vcf,
                   csv_file,
                   opt.refset,
                   opt.model,
                   gender,
                   opt.sample_id,
                   settings.vcf_template
                  ))
        print (action)
        os.system(action)

        """Make IGV session xml """ 
        if opt.model == "HC": ## Need to re-write to Ninja for proper implementation? Template + add model 'block' to xml 
            action = ("python {0} --bam {1} -o {2} --sampleid {3} --template {4} --refdate {5} --runid {6}".format(
                       settings.igv_xml,
                       opt.input_bam,
                       output_folder,
                       opt.sample_id,
                       settings.igv_xml,
                       opt.refset,
                       opt.run_id
                      ))
            print(action)
            os.system(action)
    else:
        sys.exit("Choose either make_ref or make_call")
