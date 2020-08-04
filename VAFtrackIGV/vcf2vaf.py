# Calculate Variant Allele Frequency
# Made with Python 3.7
# Required args:
#   1: path to input VCF file,
#   2: min depth
# Output BED file gets created in current work directory

import os
from argparse import ArgumentParser

def main():
    parser = ArgumentParser(description='Calculate Variant Allele Frequency')
    parser.add_argument("-i", dest="filepath", required=True,
                        help="Path to input VCF file",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-d", dest="min_depth", type=int, default=1,
                        help='Minimum depth required per variant (default: 1)')
    args = parser.parse_args()

    process_file(args.filepath, args.min_depth)

def process_file(filepath, min_depth):
    outname = os.path.splitext(os.path.basename(filepath))[0] + ".bed"
    outfile = open(outname, "w")
    outfile.write("chromosome\tstart\tstop\tVAF\n")

    with open(filepath) as file:
        header = {}
        for num, line in enumerate(file):
            if line.startswith("##"):
                continue
            elif line.startswith("#"):
                header_list = line.lstrip("#").rstrip().split("\t")
                header = {k: v for v, k in enumerate(header_list)}
            else:
                line = line.split("\t")
                info = line[-1].strip("\n").split(":")
                format = line[header["FORMAT"]].split(":")
                ref = line[header["REF"]]
                alt = line[header["ALT"]]
                # only continue with lines that are for sure SNVs,
                # filtering can be improved based on needs
                if len(ref) == 1 and len(alt) == 1:
                    if "AD" in format:
                        # find which column contains allelic depth, and retrieve
                        index = format.index("AD")
                        ad = info[index].split(",")
                        ad_ref = float(ad[0])
                        ad_alt = float(ad[1])

                        # filter on given min depth
                        if (ad_alt + ad_ref) >= min_depth:
                            vaf = (ad_alt / (ad_alt + ad_ref)) * 100
                            outfile.write("{}\t{}\t{}\t{:.2f}\n"
                                          .format(line[header["CHROM"]],
                                            int(line[header["POS"]]) - 1,  # -1 cause bed is 0-based?
                                            line[int(header["POS"])],  # length is 1
                                            vaf))
    outfile.close()
    print(f"\nDone: {outname}")


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg


if __name__ == '__main__':
    main()

