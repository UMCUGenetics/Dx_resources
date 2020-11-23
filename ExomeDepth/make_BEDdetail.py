#! /usr/bin/env python3
import os
import subprocess
import argparse
from string import Template
import statistics
import vcf
import settings

def make_merge_dic(merge_samples):
    merge_dic = {}
    with open(merge_samples, 'r') as merge_file:
        for line in merge_file:
            splitline = line.split()
            if splitline[0] not in merge_dic:
                merge_dic[splitline[0]] = []
            merge_dic[splitline[0]].append(splitline[1])

    return merge_dic

def slice_vcf(args, merge_dic):
    vcf_files = subprocess.getoutput("find {input} -type f -iname \"*vcf\"".format(input=args.inputfolder)).split()
    all_event_file = open("{output}/{outputfile}_all_events.txt".format(
        output=args.outputfolder,
        outputfile=args.outputfile
        ),"w")
    all_event_file.write("sampleID\tchromosome\tstart\tstop\tgender\trefset\tcalltype\tntargets\tsvlen\tratio\tccn\tbf\tcorrelation\tdeldupratio\ttotalcalls\trefsamples\n")
    excluded_samples_file = open("{output}/{outputfile}_excluded_samples.txt".format(
        output=args.outputfolder,
        outputfile=args.outputfile
        ),"w")

    event_dic = {}
    childs = 0
    parents = 0 
    for vcf_file in vcf_files:
        exclude = False
        sampleid = vcf_file.split("/")[-1].split("bam")[0].split("_")[2].rstrip(".")  
        runid = "_".join(vcf_file.split("/")[-1].split("bam")[1].split("_")[1:-2])

        if sampleid in merge_dic:  # Check if sample in specific run is merge sample. If so, exclude
            for run in merge_dic[sampleid]:
                if run == runid:
                    excluded_samples_file.write("Sample {sampleid} run {runid} is excluded being merge sample\n".format(
                        sampleid=sampleid,
                        runid=runid
                        ))
                    exclude = True

        if "giab" in sampleid.lower() or "control" in sampleid.lower(): # Remove GIAB and Control samples as these should not be included in the results
            excluded_samples_file.write("Sample {sampleid} run {runid} is excluded being GIAB or Control sample\n".format(
                sampleid=sampleid,
                runid=runid
                ))
            exclude = True

        if exclude == False:
            with open(vcf_file, 'r') as vcf_input_file:
                vcf_reader = vcf.Reader(open(vcf_file, 'r'))
                EDreference = vcf_reader.metadata['EDreference'][0]
                gender = EDreference.split("_")[1]
                refset = EDreference.split("_")[2]

                if len(vcf_reader.samples) > 1:  # multisample VCF
                    excluded_samples_file.write("VCF {vcf_file} from run {runid} was excluded because this was a multisample VCF\n".format(
                        vcf_file=vcf_file,
                        runid=runid
                    ))
                    break
                sampleid_vcf = vcf_reader.samples[0]

                if sampleid != sampleid_vcf: # Check if sampleID in VCF is same as sampleID in VCF. If not: report and ignore
                    excluded_samples_file.write("Sample {sampleid} run {runid} is excluded as sampleID of VCF file ({sampleid}) file is not the same as sampleID within VCF ({sampleid_vcf})\n".format(
                    sampleid=sampleid,
                    runid=runid,
                    sampleid_vcf=sampleid_vcf
                    ))
                    continue

                for record in vcf_reader:
                    """ General fields """
                    chrom = record.CHROM
                    start = record.POS

                    """ Info Field """
                    calltype = record.INFO['SVTYPE'] 
                    ntargets = record.INFO['NTARGETS']
                    svlen = record.INFO['SVLEN']
                    stop = record.INFO['END']
 
                    """ Format fields """
                    correl = float(record.genotype(sampleid_vcf)['CR'])
                    deldupratio = float(record.genotype(sampleid_vcf)['PD'])
                    totalcalls = int(record.genotype(sampleid_vcf)['TC'])
                    refsamples = record.genotype(sampleid_vcf)['RS']
                    ratio = record.genotype(sampleid_vcf)['RT']
                    ccn = record.genotype(sampleid_vcf)['CCN']
                    bf = record.genotype(sampleid_vcf)['BF']

                    if (totalcalls < args.totalcallsqc_min or
                        totalcalls > args.totalcallsqc_max or
                        correl < args.correlqc or
                        deldupratio < args.deldupratioqc_min or
                        deldupratio > args.deldupratioqc_max):
                            excluded_samples_file.write("Sample {sampleid} run {runid} is excluded because not all QC are above threshold\n".format(
                                sampleid=sampleid,
                                runid=runid
                                ))
                            break

                    all_event_file.write("{sampleid}\t{chrom}\t{start}\t{stop}\t{gender}\t{refset}\t{calltype}\t{ntargets}\t{svlen}\t{ratio}\t{ccn}\t{bf}\t{correl}\t{deldupratio}\t{totalcalls}\t{refsamples}\n".format(
                        sampleid=sampleid,
                        chrom=chrom,
                        start=start,
                        stop=stop,
                        gender=gender,
                        refset=refset,
                        calltype=calltype,
                        ntargets=ntargets,
                        svlen=svlen,
                        ratio=ratio,
                        ccn=ccn,
                        bf=bf,
                        correl=correl,
                        deldupratio=deldupratio,
                        totalcalls=totalcalls,
                        refsamples=refsamples
                        ))
 

                    event = "{chrom}_{start}_{stop}_{calltype}_{ntargets}".format(
                        chrom=chrom,
                        start=start,
                        stop=stop,
                        calltype=calltype,
                        ntargets=ntargets
                        )

                    if event not in event_dic:
                        event_dic[event] = {"parent":{"count":0, "bf":[],"ratio":[],"correlation":[],"deldupratio":[],"totalcalls":[], "gender":[]}, 
                            "child":{"count":0, "bf":[],"ratio":[],"correlation":[],"deldupratio":[],"totalcalls":[], "gender":[]}}

                    """ Determine Child or Parent status based on sampleid. There is no other option at the moment """
                    if "CM" in sampleid or "CF" in sampleid or "CO" in sampleid:
                        sampletype = "child"
                        childs +=1
                    elif "PM" in sampleid or "PF" in sampleid or "PO" in sampleid:
                        sampletype = "parent"
                        parents += 1
               
                    event_dic[event][sampletype]["count"] += 1
                    event_dic[event][sampletype]["bf"].append(bf)
                    event_dic[event][sampletype]["ratio"].append(ratio)
                    event_dic[event][sampletype]["correlation"].append(correl)
                    event_dic[event][sampletype]["deldupratio"].append(deldupratio)
                    event_dic[event][sampletype]["totalcalls"].append(totalcalls)
                    event_dic[event][sampletype]["gender"].append(gender)  #  Is not being used in the BED file output at the moment.
    all_event_file.close()
    excluded_samples_file.close()
    return event_dic, childs, parents

def make_bed_detail(args, event_dic, childs, parents):
    event_file = open("{outputfolder}/{outputfile}_UCSC.bed".format(outputfolder=args.outputfolder, outputfile=args.outputfile),"w")
    event_file_igv = open("{outputfolder}/{outputfile}_IGV.bed".format(outputfolder=args.outputfolder, outputfile=args.outputfile),"w")

    """ Write header in BED file """
    event_file.write("track name=\"HC_WES_CNV\" type=\"bedDetail\" description=\"CNVs called by Exomedepth using HC callset. #Child={childs} #Parents={parents} \" visibility=3 itemRgb=\"On\"\n".format(childs=childs, parents=parents)) 
    event_file_igv.write("track name=\"HC_WES_CNV\" type=\"bed\" description=\"CNVs called by Exomedepth using HC callset. #Child={childs} #Parents={parents} \" visibility=3 itemRgb=\"On\"\n".format(childs=childs, parents=parents))
    total_event_list = []
    for item in event_dic:
        chrom, start, stop, calltype, ntargets = item.split("_")
        """ Make start 0-based for BED file """
        start = int(start) - 1
        event_list = []
        event_list.append(chrom)
        event_list.append(start)  # Start position event
        event_list.append(stop)  # End position event
        total_count = int(event_dic[item]["parent"]["count"]) + int(event_dic[item]["child"]["count"])
        total_ratio = event_dic[item]["parent"]["ratio"] + event_dic[item]["child"]["ratio"]
        total_bf = event_dic[item]["parent"]["bf"] + event_dic[item]["child"]["bf"]
        total_ratio = [float(i) for i in total_ratio]
        total_bf = [float(i) for i in total_bf]
        median_ratio = "%.2f" % (float(statistics.median(total_ratio)))
        median_bf = "%.0f" % (float(statistics.median(total_bf)))
        summary = "{total_count}x:RT={median_ratio}:BF={median_bf}:NT={ntargets}".format(
            total_count=total_count,
            median_ratio=median_ratio,
            median_bf=median_bf,
            ntargets=ntargets
            )
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
            bfs_child = ' '.join(map(str, [elem for elem in bfs_child_list])) # Then make string to write in output
            ratios_child = ' '.join(map(str, [elem for elem in ratios_child_list]))  # Then make string to write in output

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
            bfs_parent = ' '.join(map(str, [elem for elem in bfs_parent_list])) # Then make string to write in output
            ratios_parent = ' '.join(map(str, [elem for elem in ratios_parent_list]))  # Then make string to write in output

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
    total_event_list.sort(key=lambda x:(settings.chromosome_order[x[0]], int(x[1])))
      
    for item in total_event_list:
        item[0] = "chr{original}".format(original=item[0])
        joined_item = "\t".join([str(i) for i in item])
        event_file.write("{joined}\n".format(joined=joined_item))
      
        joined_item_igv = "\t".join([str(i) for i in item[0:-1]])
        event_file_igv.write("{joined}\n".format(joined=joined_item_igv))

    event_file.close()
    event_file_igv.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Path to folder including VCF files')
    parser.add_argument('outputfolder', help='Path to output folder')
    parser.add_argument('outputfile', help='output prefix filename')
    parser.add_argument('merge_samples', help='Path to file including all merge samples (tab delimited file with column 1 = sampleID, 2 = run/projectID')
    parser.add_argument('--totalcallsqc_min', default=settings.number_calls[0], type=int, help='Threshold for minimum allowed CNV calls (default = 35)')
    parser.add_argument('--totalcallsqc_max', default=settings.number_calls[1], type=int, help='Threshold for maximum allowed CNV calls (default = 200)')
    parser.add_argument('--deldupratioqc_min', default=settings.del_dup_ratio[0], type=float, help='Threshold for minimum allowed deldupratio (default = 15)')
    parser.add_argument('--deldupratioqc_max', default=settings.del_dup_ratio[1], type=float, help='Threshold for maximum allowed deldupratio (default = 85)')
    parser.add_argument('--correlqc', default=settings.correlation, type=float, help='Threshold for minimum allowed correlation score (default = 0.98)')
    args = parser.parse_args()

    if not os.path.isdir(args.outputfolder):
        os.system("mkdir -p {outputfolder}".format(outputfolder=args.outputfolder))
    merge_dic = make_merge_dic(args.merge_samples)
    event_dic, childs, parents = slice_vcf(args, merge_dic)
    make_bed_detail(args, event_dic, childs, parents)
