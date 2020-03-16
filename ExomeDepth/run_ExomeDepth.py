#! /usr/bin/env python
import os
import sys
import pysam
import commands
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
    group.add_option("-i", dest = "input_folder",
                     metavar = "[PATH]", help = "path to input folder folder")
    group.add_option("--ib", dest = "input_bam", 
                     metavar = "[PATH]", help = "path to input bam file")
    group.add_option("-o", dest = "output_folder",
                     default = "./", metavar = "[PATH]",
                     help = "path to output folder [default = ./]")
    group.add_option("-g", dest = "gender_file",
                     metavar = "[PATH]",
                     help = "full path to tab delimited txt file with bam \
                    + gender (male/female) [default = off]")
    group.add_option("-p", dest = "prefix", metavar = "[PATH]",
                     help = "referenceset naming [default = input_folder basename]")
    group.add_option("-m", dest = "mail", metavar = "[STRING]",
                     help = "email adress of submitter [default = None")
    parser.add_option_group(group)
    (opt, args) = parser.parse_args()

if opt.input_folder and opt.input_bam:
    sys.exit("choose either folder or bam file")

if opt.input_folder:
    wkdir = opt.input_folder.rstrip("/")
elif opt.input_bam:
    wkdir = "/".join(opt.input_bam.split("/")[0-1]).rstrip("/")
else:
    wkdir = commands.getoutput("pwd").rstrip("/")

outdir = opt.output_folder.rstrip("/")
if not opt.mail and not opt.input_bam:
    sys.exit("provide email")

if not os.path.isdir(outdir):
    os.system("mkdir -p " + str(outdir))

if opt.prefix:
    prefix = opt.prefix
else:
    prefix = str(wkdir).split("/")[-1]

if opt.gender_file:
    gender_dic = {}
    gender_file = open(str(opt.gender_file), "r")
    for line in gender_file:
        splitline = line.split()
        if splitline[0] not in gender_dic:
            gender_dic[splitline[0]] = splitline[1]
    gender_file.close()

analysis = settings.analysis
ed_r = settings.call_cnv_r
refgenome = settings.reference_genome
qsub_call = str(settings.qsub_call) + str(opt.mail)
qsub_ref = str(settings.qsub_ref) + str(opt.mail)
prob = settings.probability
gender = settings.gender

"""Log all settings in setting.log file"""
log_dir = str(outdir)+"/logs"
os.system("mkdir -p "+str(log_dir))
if opt.input_folder:
    write_file = open(log_dir+"/settings.log", "w")
elif opt.input_bam:
    log_file = "{0}/{1}_settings.log".format(str(log_dir), str(opt.input_bam.split("/")[-1]))
    write_file = open(log_file, "w") 

(options, args) = parser.parse_args()
for item in vars(options):
    write_file.write("{0}\t{1}\n".format(str(item), str(vars(options)[item])))

for item in dir(settings):
    if "__" not in item:
        write_file.write("{0}\tt{1}\n".format(item, str(repr(eval("settings.%s" % item)))))
write_file.close()

if opt.make_ref and not opt.make_call:
    """Make new reference set."""
    if opt.input_bam:
        sys.exit("please provide BAM folder")

    bams =  commands.getoutput("find -L {0} -iname \"*realigned.bam\" ".format(wkdir)).split()
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
            folder = "{0}/{1}_{2}_{3}".format(outdir, target, item, str(wkdir).split("/")[-1])  
            output_id = "{0}_{1}_{2}.EDref".format(target, item, prefix)
            os.system("mkdir -p {0}".format(folder))
            for bam in ref_gender_dic[item]:
                os.system("ln -sd " + str(bam) + "* " + str(folder))
            write_file = open(str(folder) + "/make_ref.sh", "w")
            write_file.write(str(qsub_ref) + "\n")
            write_file.write("module load {0}\n".format(settings.r_version))
            write_file.write("cd {0}\n".format(folder))
            write_file.write("Rscript {0} {1}/ {1}/{2} {3} {4} {5}\n".format(
                             settings.create_refset_r,
                             folder,
                             output_id,
                             analysis[target]["target_bed"],
                             settings.reference_genome,
                             analysis[target]["exon_bed"]
                             ))
            write_file.close()
            os.chdir(str(folder))
            os.system("qsub {0}/make_ref.sh".format(folder))

elif opt.make_call and not opt.make_ref:  # Call CNV from BAMs
    """Make CNV call on sample(s)."""

    if opt.input_bam:
        bams = [str(opt.input_bam)]
    else:
        bams = commands.getoutput("find -L {0} -iname \"*realigned.bam\"".format(wkdir)).split()
    print("Number of BAM files detected = {0}".format(len(bams)))
    if not bams:
        sys.exit("no bams detected")

    for bam in bams:
        """Get gender from chrY read count ratio."""
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
                print("Sample {0} has a unknown gender and will not be analysed".format(bam.split("/")[-1]))
                continue

        for item in analysis:
            print("Submitting {0} jobs".format(item))
            sampleid = bam.split("/")[-1].split("_")[0]
            outfolder = "{0}/{1}/{1}_{2}".format(outdir, item, bam.split("/")[-1])
            os.system("mkdir -p {0}".format(outfolder))
            if opt.input_bam:  #Single BAM processing in serial
                os.chdir(outfolder)
                os.system("module load {0} && Rscript {1} {2} {3} {4} {5} {6} {7}".format(
                          settings.r_version,
                          ed_r,
                          analysis[item]["refset"][gender],
                          analysis[item]["target_bed"],
                          refgenome,
                          analysis[item]["exon_bed"],
                          str(prob[str(item)]),
                          bam
                         ))
                os.system("rename {0} {1}_{2}_{3} *".format(
                          bam.split("/")[-1],
                          item,
                          settings.refset,
                          bam.split("/")[-1]
                          ))
                os.system("python {0} -i {1} -t {2} -m {3} --id={4}".format(
                          settings.csv2vcf,
                          outfolder,
                          settings.template,
                          analysis[item]["refset"][gender],
                          sampleid
                         ))
                os.system("cp {0}/*vcf {1}/{2}".format(
                          outfolder,
                          outdir,
                          item
                         ))
            else:  # multiple BAM processing for all BAMs in folder (SGE)
                write_file = open("{0}/{1}_{2}_{3}.sh".format(
                                outfolder,
                                item,
                                gender,
                                bam.split("/")[-1][:-4]		#[:-4] easier for renaming before csv2vcf
                                ), "w")  
                write_file.write(str(qsub_call) + "\n")
                write_file.write("module load {0}\n".format(settings.r_version))
                write_file.write("cd {0}\n".format(outfolder))
                write_file.write("Rscript {0} {1} {2} {3} {4} {5} {6}\n".format(
                                 ed_r,
                                 analysis[item]["refset"][gender],
                                 analysis[item]["target_bed"],
                                 refgenome,
                                 analysis[item]["exon_bed"],
                                 str(prob[str(item)]),
                                 bam
                                 ))
                write_file.close()
                """Submit jobs."""
                os.chdir(outfolder)
                command = "qsub {0}/{1}_{2}_{3}.sh ".format(
                         outfolder,
                         item,
                         gender,
                         bam.split("/")[-1][:-4]
                         ) 
                job_id = commands.getoutput(command).split()[2]
                """Make VCF per sample here with ed_csv_to_vcf.py."""
                write_file = open("{0}/{1}_{2}_csv_2_vcf.sh".format(
                                  outfolder,
                                  item,
                                  gender
                                 ),"w")
                write_file.write(str(qsub_call) + "\n")
                write_file.write("#$ -hold_jid {0}\n".format(job_id))
                write_file.write("cd {0}\n".format(outfolder))
                write_file.write("rename {0} {1}_{2}_{3} *\n".format(
                                 bam.split("/")[-1],
                                 item,
                                 settings.refset,
                                 bam.split("/")[-1]
                                 ))
                write_file.write("python {0} -i {1} -t {2} -m {3} --id={4}\n".format(
                                 settings.csv2vcf,
                                 outfolder,
                                 settings.template,
                                 analysis[item]["refset"][gender],
                                 sampleid
                                 ))
                write_file.write("cp {0}/*vcf {1}/{2} \n".format(
                                 outfolder,
                                 outdir,
                                 item
                                 ))
                write_file.close()
                os.system("qsub {0}/{1}_{2}_csv_2_vcf.sh".format(outfolder, item, gender))
        """Touch done file is loop is completed"""
        os.system("touch {0}/{1}.done".format(log_dir, bam.split("/")[-1]))
else:
    sys.exit("Choose either make_ref or make_call")
