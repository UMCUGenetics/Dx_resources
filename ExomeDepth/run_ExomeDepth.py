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
    if (read.mapping_quality >= 20 and
            read.reference_end and read.reference_start is not None):
        return True
    else:
        return False


def get_gender(bam):
    """Determine chrY ratio based on read count in bam (excl PAR)."""
    workfile = pysam.AlignmentFile(bam, "rb")
    locus = 'Y:2649520-59034050'
    yreads = float(sum([valid_read(read)
                        for read in workfile.fetch(region=locus)
                        ]
                       )
                   )
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
                     action='store_true', help="Make new reference set")
    group.add_option("-c", dest="make_call",
                     action='store_true', help="Call CNV from BAMs")
    group.add_option("-i", dest="input_folder",
                     metavar="[PATH]", help="path to input folder folder")
    group.add_option("--ib", dest="input_bam", 
                     metavar="[PATH]", help="path to input bam file")
    group.add_option("-o", dest="output_folder",
                     default="./", metavar="[PATH]",
                     help="path to output folder [default = ./]")
    group.add_option("-g", dest="gender_file",
                     metavar="[PATH]",
                     help="full path to tab delimited txt file with bam \
                    + gender (male/female) [default = off]")
    group.add_option("-p", dest="prefix", metavar="[PATH]",
                     help="referenceset naming [default = input_folder basename]")
    group.add_option("-m", dest="mail", metavar="[STRING]",
                     help="email adress of submitter [default = None")
    parser.add_option_group(group)
    (opt, args) = parser.parse_args()

if opt.input_folder and opt.input_bam:
    sys.exit("choose either folder or bam file")

if opt.input_folder:
    wkdir = opt.input_folder.rstrip("/")
elif opt.input_bam:
    wkdir="/".join(opt.input_bam.split("/")[0-1]).rstrip("/")
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
    g_dic = {}
    g_file = open(str(opt.gender_file), "r")
    for line in g_file:
        splitline = line.split()
        if splitline[0] not in g_dic:
            g_dic[splitline[0]] = splitline[1]

analysis = settings.analysis
ed_r = settings.callcnv_r
refgenome = settings.reference_genome
qsub_call = str(settings.qsub_call) + str(opt.mail)
qsub_ref = str(settings.qsub_ref) + str(opt.mail)
prob = settings.probability
gender = settings.gender


"""Log all settings in setting.log file"""
log_dir=str(outdir)+"/logs"
os.system("mkdir -p "+str(log_dir))
if opt.input_folder:
    write_file=open(log_dir+"/settings.log","w")
elif opt.input_bam:
     write_file=open(log_dir+"/"+str(opt.input_bam.split("/")[-1])+"_settings.log","w")
(options, args) = parser.parse_args()
for item in vars(options):
    write_file.write(str(item)
                     + "\t"
                     + str(vars(options)[item])
                     + "\n"
                     )
for item in dir(settings):
    if "__" not in item:
        write_file.write(str(item)
                         + "\t"
                         + str(repr(eval("settings.%s" % item)))
                         + "\n"
                         )
write_file.close()

if opt.make_ref and not opt.make_call:
    """Make new reference set."""

    if opt.input_bam:
        sys.exit("please provide BAM folder")

    bams = commands.getoutput("find -L "
                              + str(wkdir)
                              + " -iname \"*realigned.bam\" "
                              ).split()
    print "Number of BAM files detected = " + str(len(bams))

    """Get gender from chrY read count ratio."""
    dic = {}
    for item in gender:
        if str(item) not in dic:
            dic[str(item)] = []

    for bam in bams:
        gender = get_gender(bam)
        if gender is not "unknown":
            dic[get_gender(bam)] += [bam]
        else:
            print ("Sample "
                   + str(bam)
                   + "  has unknown gender and is removed from analysis"
                   )

    """Make folder per gender + analysis, and soflink BAMs in these folders."""
    for target in analysis:
        for item in dic:
            folder = (str(outdir) + "/" + str(target) + "_"+str(item) + "_"
                      + str(wkdir).split("/")[-1]
                      )
            output_id = (str(target) + "_" + str(item) + "_"
                         + str(prefix) + ".EDref"
                         )
            os.system("mkdir -p " + str(folder))
            for bam in dic[item]:
                os.system("ln -sd " + str(bam) + "* " + str(folder))
            write_file = open(str(folder) + "/make_ref.sh", "w")
            write_file.write(str(qsub_ref) + "\n")
            write_file.write("module load " + str(settings.Rversion)+"\n")
            write_file.write("cd " + str(folder) + "\n")
            write_file.write("Rscript " + str(settings.createrefset_r) + " "
                             + str(folder) + "/ " + str(folder) + "/"
                             + str(output_id) + " "
                             + str(analysis[target]["target_bed"]) + " "
                             + str(settings.reference_genome) + " "
                             + str(analysis[target]["exon_bed"]) + "\n"
                             )
            write_file.close()
            os.chdir(str(folder))
            os.system("qsub " + str(folder) + "/make_ref.sh")

elif opt.make_call and not opt.make_ref:  # Call CNV from BAMs
    """Make CNV call on sample(s)."""

    if opt.input_bam:
        bams=[str(opt.input_bam)]
    else:
        bams= commands.getoutput("find -L " + str(wkdir) 
                                 + " -iname \"*realigned.bam\" "
                                 ).split()

    print "Number of BAM files detected = " + str(len(bams))
    if not bams:
        sys.exit("no bams detected")

    for bam in bams:
        """Get gender from chrY read count ratio."""
        gender = get_gender(bam)
        if opt.gender_file:  # overrule gender as given in gender_file
            for item in g_dic:
                if str(item) in str(bam):
                    gender = g_dic[item]
        if gender == "unknown":  # try to determine gender on BAM ID annotation
            if re.search('[C|P]M', bam.split("/")[-1]):
                gender = "male"
            elif re.search('[C|P]F', bam.split("/")[-1]):
                gender = "female"
            else:
                print ("Sample " + str(bam.split("/")[-1])
                       + " has a unknown gender and will not be analysed"
                       )
                continue

        for item in analysis:
            print "submitting " + str(item) + " jobs"
            sampleid = bam.split("/")[-1].split("_")[0]
            outfolder = (str(outdir) + "/"+ str(item)+"/" + str(item) + "_"
                         + str(bam.split("/")[-1])
                         )
            os.system("mkdir -p " + str(outfolder))
            if opt.input_bam:  #Single BAM processing in serial
                os.chdir(outfolder)
                os.system("module load " + str(settings.Rversion) + "&& Rscript " 
                          + str(ed_r) + " " + str(analysis[item]["refset"][gender]) 
                          + " " + str(analysis[item]["target_bed"]) + " " 
                          + str(refgenome) + " " + str(analysis[item]["exon_bed"]) 
                          + " " + str(prob[str(item)]) + " "+str(bam)
                          )
                os.system("rename " + str(bam.split("/")[-1]) + " " 
                          + str(item)+"_" + str(settings.refset)
                          + "_" + str(bam.split("/")[-1]) + " *"
                          )
                os.system("python " + str(settings.csv2vcf) + " -i " + str(outfolder) 
                          + " -t " + str(settings.template) + " -m " 
                          + str(analysis[item]["refset"][gender]) 
                          + " --id "+ str(sampleid)
                          )
                os.system("cp " + str(outfolder) + "/*vcf " + str(outdir) + "/" + str(item))

            else:  # multiple BAM processing for all BAMs in folder (SGE)
                write_file = open(str(outfolder) + "/" + str(item)
                                  + "_" + str(gender) + "_"
                                  + bam.split("/")[-1][:-4] + ".sh", "w"	
                                  )   #[:-4] easier for renaming before csv2vcf					
                write_file.write(str(qsub_call) + "\n")
                write_file.write("module load " + str(settings.Rversion) + "\n")
                write_file.write("cd " + str(outfolder) + "\n")
                write_file.write("Rscript " + str(ed_r) + " "
                                 + str(analysis[item]["refset"][gender])
                                 + " " + str(analysis[item]["target_bed"])
                                 + " " + str(refgenome) + " "
                                 + str(analysis[item]["exon_bed"])
                                 + " " + str(prob[str(item)]) 
                                 + " " + str(bam) + "\n"
                                 )
                write_file.close()
                """Submit jobs."""
                os.chdir(outfolder)
                command = ("qsub " + str(outfolder) + "/"
                           + str(item) + "_" + str(gender)
                           + "_" + bam.split("/")[-1][:-4] + ".sh"
                           )
                job_id = commands.getoutput(command).split()[2]
                """Make VCF per sample here with ed_csv_to_vcf.py."""

                write_file = open(outfolder + "/"
                                  + str(item) + "_" + str(gender)
                                  + "_" + "csv_2_vcf.sh", "w"
                                  )
                write_file.write(str(qsub_call) + "\n")
                write_file.write("#$ -hold_jid " + str(job_id) + "\n")
                write_file.write("cd " + str(outfolder)+ "\n")
                write_file.write("rename " + str(bam.split("/")[-1]) + " " 
                                 + str(item)+"_" + str(settings.refset) 
                                 + "_" + str(bam.split("/")[-1]) + " *\n"
                                 )
                write_file.write("python " + str(settings.csv2vcf) + " -i "
                                 + str(outfolder) + " -t " + str(settings.template)
                                 + " -m " + str(analysis[item]["refset"][gender])
                                 + " --id " + str(sampleid) + "\n"
                                 )
                write_file.write("cp " + str(outfolder) + "/*vcf " + str(outdir) + "/" + str(item) + "\n")
                write_file.close()
                os.system("qsub " + outfolder + "/" + str(item) + "_"
                          + str(gender) + "_" + "csv_2_vcf.sh"
                          )
        """Touch done file is loop is completed"""
        os.system("touch "+str(log_dir) + "/"+str(bam.split("/")[-1])+".done")
else:
    sys.exit("Choose either make_ref or make_call")
