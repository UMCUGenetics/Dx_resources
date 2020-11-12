#! /usr/bin/env python3
import os
import subprocess
import argparse
import statistics
from string import Template
import settings

def make_merge_dic(merge_samples):
    merge_dic = {}
    with open(merge_samples, 'r') as merge_file:
        for line in merge_file:
            splitline = line.split()
            if splitline[0] not in merge_dic:
                merge_dic[splitline[0]] = [splitline[1]]
            else:
                merge_dic[splitline[0]].append(splitline[1])
        return merge_dic

def slice_vcf (args, merge_dic):
    vcf_files = subprocess.getoutput("find {} -type f -iname \"*vcf\"".format(args.inputfolder)).split()
    event_file = open("{}/{}".format(args.outputfolder,"All_events.txt"),"w")
    event_file.write("sampleid\tchromsome\tstart\tstop\tgender\tcalltype\tntargets\tsvlen\tratio\tccn\tbf\tcorrelation\tdeldupratio\ttotalcalls\rrefsamples\n")
    event_dic = {} 
    for vcf in vcf_files:
        exclude = False
        sampleid = vcf.split("/")[-1].split("bam")[0].split("_")[2].rstrip(".")  
        runid = "_".join(vcf.split("/")[-1].split("bam")[1].split("_")[1:-2])

        if sampleid in merge_dic:
            for run in merge_dic[sampleid]:
                if run == runid:
                    print("Sample {} run {} is excluded being merge sample".format(sampleid,runid))
                    exclude = True

        if "giab" in sampleid.lower() or "control" in sampleid.lower():
            print("Sample {} run {} is excluded being GIAB or Control sample".format(sampleid,runid))
            exclude = True

        if exclude == False:
            with open(vcf, 'r') as vcf_lines:
                """ Extract relevant fields from lines """
                for line in vcf_lines:
                    splitline = line.split()
                    if "##" in line:
                        if "EDreference" in line: 
                            gender = line.split("=")[1].split("_")[1].rstrip()
                            refset = line.split("=")[1].split("_")[2].rstrip()
                    elif "#" in line:
                        sampleid_vcf = str(splitline[9].rstrip())
                        if sampleid != sampleid_vcf: # Check if sampleID in VCF is same as sampleID in VCF. If not: report and ignore
                            print ("Sample {} run {} is excluded as sampleID of VCF file ({}) file is not the same as sampleID within VCF ({})".format(sampleid,runid,sampleid,sampleid_vcf)) 
                            break
                    else:   
                        infofields = splitline[7].split(";")
                        formatfields = splitline[9].split(":")

                        """ Sample stats """
                        correl = formatfields[4]
                        deldupratio = formatfields[8]
                        totalcalls = formatfields[9]
                        refsamples = formatfields[5] 
                        if int(totalcalls) < int(args.totalcallsqc_min) or int(totalcalls) > int(args.totalcallsqc_max) or float(correl) < float(args.correlqc) or float(deldupratio) < float(args.deldupratioqc_min) or float(deldupratio) > float(args.deldupratioqc_max):
                            print("Sample {} run {} is excluded because not all QC are above threshold".format(sampleid,runid))
                            break

                        """ Event stats """
                        chrom = splitline[0]
                        start = splitline[1]
                        stop = infofields[2].split("=")[1] 
                        calltype = infofields[0].split("=")[1]
                        ntargets = infofields[1].split("=")[1]
                        svlen = infofields[3].split("=")[1]
                        ratio = formatfields[3]
                        ccn = formatfields[1]
                        bf = formatfields[2]
                        event_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            sampleid, chrom, start, stop, gender, refset, calltype,ntargets,svlen,ratio,ccn,bf,correl,deldupratio,totalcalls,refsamples
                        ))
                        event = "{}_{}_{}_{}_{}".format(chrom, start, stop, calltype,ntargets)
                        if event not in event_dic:
                            event_dic[event] = {"parent":{"count":0, "bf":[],"ratio":[],"correlation":[],"deldupratio":[],"totalcalls":[], "gender":[]}, \
                                            "child":{"count":0, "bf":[],"ratio":[],"correlation":[],"deldupratio":[],"totalcalls":[], "gender":[]}}

                        gender_dic = {"female":"F", "male":"M"}
                        if "CM" in sampleid or "CF" in sampleid or "CO" in sampleid:
                            sampletype = "child"
                        elif "PM" in sampleid or "PF" in sampleid or "PO" in sampleid:
                            sampletype = "parent"
                
                        event_dic[event][sampletype]["count"] += 1
                        event_dic[event][sampletype]["bf"].append(bf)
                        event_dic[event][sampletype]["ratio"].append(ratio)
                        event_dic[event][sampletype]["correlation"].append(correl)
                        event_dic[event][sampletype]["deldupratio"].append(deldupratio)
                        event_dic[event][sampletype]["totalcalls"].append(totalcalls)
                        event_dic[event][sampletype]["gender"].append(gender_dic[gender])
    event_file.close()
    return event_dic

def make_beddetail(args, event_dic):
    event_file = open("{}/{}".format(args.outputfolder,args.outputfile),"w")
    """ print header in BED file """
    event_file.write("track name=\"HC_WES_CNV\" type=\"bedDetail\" description=\"CNVs called by Exomdepth using HC callset\" visibility=3 itemRgb=\"On\"\n") 
    total_event_list = []
    for item in event_dic:
        chrom, start, stop, calltype, ntargets = item.split("_")
        """ Make start 0-based for BED file """
        start = int(start) - 1
        event_list = []

        if chrom == "X":
            event_list.append(23)
        elif chrom == "Y":
            event_list.append(24)
        else:
            event_list.append(int(chrom))
        event_list.append(start)  # Start position event
        event_list.append(stop)  # End position event
        total_count = int(event_dic[item]["parent"]["count"]) + int(event_dic[item]["child"]["count"])
        total_ratio = event_dic[item]["parent"]["ratio"] + event_dic[item]["child"]["ratio"]
        total_bf = event_dic[item]["parent"]["bf"] + event_dic[item]["child"]["bf"]
        total_ratio = [float(i) for i in total_ratio]
        total_bf = [float(i) for i in total_bf]
        median_ratio = "%.2f" % (float(statistics.median(total_ratio)))
        median_bf = "%.0f" % (float(statistics.median(total_bf)))
        summary = "{}x:RT={}:BF={}:NT={}".format(total_count,median_ratio,median_bf,ntargets)
        event_list.append(summary)  # Name

        """ Append unused but necessary columns for bed file """
        event_list.append(0)  # Score
        event_list.append(".")  # Strand
        event_list.append(start)  # thickStart
        event_list.append(stop)  # thickEnd

        """ Append color code based on DEL (red) or DUP (blue) """
        if float(median_ratio) >= 1:  # Assumed DUP
            color = "0,0,255"  # itemRgb
        else:  # Assumed DEL
            color = "255,0,0"  # itemRgb
        event_list.append(color)

        """ Append blockCount, blockSizes, blockStarts. Note that these are not used for ntargets as start/stop for exons is unknown"""
        event_list.append(1)  # blockCount
        event_list.append(int(stop)-start)  # blockSizes
        event_list.append(0)  # blockStarts.

        """ Append custom annotation field in BEDdetail format. 2 additional colums are allowed! """
        """ custom1 =  sample specific (Correlation and #Calls, split for child and parent). Custom2 = event specific (Ratio and BF) """
        correlation_child = [float(i) for i in event_dic[item]["child"]["correlation"]]
        correlation_parent = [float(i) for i in event_dic[item]["parent"]["correlation"]]
        deldupratio_child = [float(i) for i in event_dic[item]["child"]["deldupratio"]]
        deldupratio_parent = [float(i) for i in event_dic[item]["parent"]["deldupratio"]]
        totalcalls_child = [float(i) for i in event_dic[item]["child"]["totalcalls"]]
        totalcalls_parent = [float(i) for i in event_dic[item]["parent"]["totalcalls"]]
        bf_child = [float(i) for i in event_dic[item]["child"]["bf"]]
        bf_parent = [float(i) for i in event_dic[item]["parent"]["bf"]]
        ratio_child = [float(i) for i in event_dic[item]["child"]["ratio"]]
        ratio_parent = [float(i) for i in event_dic[item]["parent"]["ratio"]]
        totalcount_child = int(event_dic[item]["child"]["count"])
        totalcount_parent = int(event_dic[item]["parent"]["count"])

        """ Calculate Child specific stats """
        average_correlation_child = "n/a"
        average_deldupratio_child = "n/a"
        average_totalcalls_child = "n/a"
        average_bf_child = "n/a"
        average_ratio_child = "n/a"
        bfs_child = "n/a"
        ratios_child = "n/a"
        stdev_correlation_child = "n/a"
        stdev_deldupratio_child = "n/a"
        stdev_totalcalls_child = "n/a"
        stdev_bf_child = "n/a"
        stdev_ratio_child = "n/a"

        if totalcount_child >= 1:
            bfs_child_list = list(map(float, ['%.1f' % elem for elem in bf_child][0:20]))   # Select first 20 BFs. No sorting yet as this would not be a random subset
            ratios_child_list = list(map(float, ['%.1f' % elem for elem in ratio_child][0:20]))  # Select first 20 (random) Ratios.  No sorting yet as this would not be a random subset
            bfs_child_list.sort() # Then sort
            ratios_child_list.sort() # Then sort
            bfs_child = ' '.join(map(str, [elem for elem in bfs_child_list])) # Then make string to print in output
            ratios_child = ' '.join(map(str, [elem for elem in ratios_child_list]))  # Then make string to print in output

        if totalcount_child == 1:
            average_correlation_child = "%.3f" % (correlation_child[0])
            average_deldupratio_child = "%.0f" % (deldupratio_child[0])
            average_totalcalls_child = "%.0f" % (totalcalls_child[0])
            average_bf_child = "%.1f" % (bf_child[0])
            average_ratio_child = "%.2f" % (ratio_child[0])
        elif totalcount_child > 1 :
            stdev_correlation_child = "%.3f" % (float(statistics.stdev(correlation_child)))
            average_correlation_child = "%.3f" % (float(statistics.mean(correlation_child)))
            stdev_deldupratio_child = "%.0f" % (float(statistics.stdev(deldupratio_child)))
            average_deldupratio_child = "%.0f" % (float(statistics.mean(deldupratio_child)))
            stdev_totalcalls_child = "%.0f" % (float(statistics.stdev(totalcalls_child)))
            average_totalcalls_child = "%.0f" % (float(statistics.mean(totalcalls_child)))
            stdev_bf_child = "%.1f" % (float(statistics.stdev(bf_child)))
            average_bf_child = "%.1f" % (float(statistics.mean(bf_child)))
            stdev_ratio_child = "%.2f" % (float(statistics.stdev(ratio_child)))
            average_ratio_child = "%.2f" % (float(statistics.mean(ratio_child)))

        """ Calculate Parent specific stats """
        average_correlation_parent = "n/a"
        average_deldupratio_parent = "n/a"
        average_totalcalls_parent = "n/a"
        average_bf_parent = "n/a"
        average_ratio_parent = "n/a"
        bfs_parent = "n/a"
        ratios_parent = "n/a"
        stdev_correlation_parent = "n/a"
        stdev_deldupratio_parent = "n/a"
        stdev_totalcalls_parent = "n/a"
        stdev_bf_parent = "n/a"
        stdev_ratio_parent = "n/a"

        if totalcount_parent >= 1:
            bfs_parent_list = list(map(float, ['%.1f' % elem for elem in bf_parent][0:20]))   # Select first 20 BFs. No sorting yet as this would not be a random subset
            ratios_parent_list = list(map(float, ['%.1f' % elem for elem in ratio_parent][0:20]))  # Select first 20 (random) Ratios.  No sorting yet as this would not be a random subset
            bfs_parent_list.sort() # Then sort
            ratios_parent_list.sort() # Then sort
            bfs_parent = ' '.join(map(str, [elem for elem in bfs_parent_list])) # Then make string to print in output
            ratios_parent = ' '.join(map(str, [elem for elem in ratios_parent_list]))  # Then make string to print in output

        if totalcount_parent == 1:
            average_correlation_parent = "%.3f" % (correlation_parent[0])
            average_deldupratio_parent = "%.0f" % (deldupratio_parent[0])
            average_totalcalls_parent = "%.0f" % (totalcalls_parent[0])
            average_bf_parent = "%.1f" % (bf_parent[0])
            average_ratio_parent = "%.2f" % (ratio_parent[0])
        elif totalcount_parent > 1 :
            stdev_correlation_parent = "%.3f" % (float(statistics.stdev(correlation_parent)))
            average_correlation_parent = "%.3f" % (float(statistics.mean(correlation_parent)))
            stdev_deldupratio_parent = "%.0f" % (float(statistics.stdev(deldupratio_parent)))
            average_deldupratio_parent = "%.0f" % (float(statistics.mean(deldupratio_parent)))
            stdev_totalcalls_parent = "%.0f" % (float(statistics.stdev(totalcalls_parent)))
            average_totalcalls_parent = "%.0f" % (float(statistics.mean(totalcalls_parent)))
            stdev_bf_parent = "%.1f" % (float(statistics.stdev(bf_parent)))
            average_bf_parent = "%.1f" % (float(statistics.mean(bf_parent)))
            stdev_ratio_parent = "%.2f" % (float(statistics.stdev(ratio_parent)))
            average_ratio_parent = "%.2f" % (float(statistics.mean(ratio_parent)))

        template_file = Template(open(settings.html).read())
        substitute_dic = {'total_child' : totalcount_child, 
            'cr_avr_child' : average_correlation_child, 
            'cr_stdev_child' : stdev_correlation_child,
            'deldup_avr_child': average_deldupratio_child, 
            'deldup_stdev_child' : stdev_deldupratio_child, 
            'numbercalls_avr_child' : average_totalcalls_child, 
            'numbercalls_stdev_child' : stdev_totalcalls_child,
            'total_parent' : totalcount_parent,
            'cr_avr_parent' : average_correlation_parent,
            'cr_stdev_parent' : stdev_correlation_parent,
            'deldup_avr_parent': average_deldupratio_parent,
            'deldup_stdev_parent' : stdev_deldupratio_parent,
            'numbercalls_avr_parent' : average_totalcalls_parent,
            'numbercalls_stdev_parent' : stdev_totalcalls_parent,
            'average_bf_child': average_bf_child, 
            'stdev_bf_child': stdev_bf_child,
            'average_ratio_child': average_ratio_child,
            'stdev_ratio_child': stdev_ratio_child,
            'average_bf_parent': average_bf_parent,
            'stdev_bf_parent': stdev_bf_parent,
            'average_ratio_parent': average_ratio_parent,
            'stdev_ratio_parent': stdev_ratio_parent,
            'bfs_child' : bfs_child,
            'ratios_child': ratios_child,
            'bfs_parent' : bfs_parent,
            'ratios_parent': ratios_parent
            }
        new_file = template_file.substitute(substitute_dic)       
        event_list.append("".join(new_file.split("\n")))
        total_event_list.append(event_list)
    total_event_list.sort(key=lambda x:(x[0], int(x[1])))  # Sort list on chromosome, start position

    for item in total_event_list:
        if item[0] == 23:
            item[0] = "chrX"
        elif item[0] == 24:
            item[0] = "chrY"
        else:
            item[0] = "chr{}".format(item[0])
        joined_item = "\t".join([str(i) for i in item])
        event_file.write("{}\n".format(joined_item))
    event_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Path to folder including VCF files')
    parser.add_argument('outputfolder', help='Path to output folder')
    parser.add_argument('outputfile', help='output filename of bed file')
    parser.add_argument('merge_samples', help='Path to file including all merge samples (tab delimited file with column 1 = sampleID, 2 = run/projectID')
    parser.add_argument('--totalcallsqc_min', default=35, help='Threshold for minimum allowed CNV calls (default = 35)')
    parser.add_argument('--totalcallsqc_max', default=200, help='Threshold for maximum allowed CNV calls (default = 200)')
    parser.add_argument('--deldupratioqc_min', default=15, help='Threshold for minimum allowed deldupratio (default = 15)')
    parser.add_argument('--deldupratioqc_max', default=85, help='Threshold for maximum allowed deldupratio (default = 85)')
    parser.add_argument('--correlqc', default=0.98, help='Threshold for minimum allowed correlation score (default = 0.98)')
    args = parser.parse_args()

    if not os.path.isdir(args.outputfolder):
        os.system("mkdir -p {}".format(args.outputfolder))
    merge_dic = make_merge_dic(args.merge_samples)
    event_dic = slice_vcf(args, merge_dic)
    make_beddetail(args, event_dic)
