#! /usr/bin/env python3

import os
import subprocess
import argparse
import random
from string import Template
import statistics
import vcf

import genologics.lims

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


def get_familystatus_clarity(sampleid, lims):
    """Get family status from Clarity LIMS."""
    if "-" in sampleid:
        sampleid = sampleid.split("-")[0]
    samples = lims.get_samples(udf={settings.monster_udf: sampleid})
    familystatus = []
    for sample in samples:
        familystatus.append(settings.family_translation[sample.udf[settings.familie_udf].lower()])
    if len(set(familystatus)) == 1:
        return familystatus[0]
    else:
        return "unknown"


def get_familystatus_sampleid(sampleid):
    if any([child_tag in sampleid for child_tag in ['CM', 'CF', 'CO']]):
        return "child"
    elif any([parent_tag in sampleid for parent_tag in ['PM', 'PF', 'PO']]):
        return "parent"
    else:
        return "unknown"


def slice_vcf(args, merge_dic, lims):
    vcf_files = subprocess.getoutput("find -L {input} -type f -iname \"*vcf\"".format(
        input=args.inputfolder
    )).split()
    all_event_file = "{output}/{outputfile}_all_events.txt".format(
        output=args.outputfolder, outputfile=args.outputfile
    )
    excluded_samples_file = "{output}/{outputfile}_excluded_samples.txt".format(
        output=args.outputfolder, outputfile=args.outputfile
    )
    with open(all_event_file, "w") as all_event_file, open(excluded_samples_file, "w") as excluded_samples_file:
        all_event_file.write(
            "sampleID\tchromosome\tstart\tstop\tgender\trefset\tcalltype\t"
            "ntargets\tsvlen\tratio\tccn\tbf\tcorrelation\tdeldupratio\ttotalcalls\trefsamples\n"
        )

        event_dic = {}
        count_familystatus = {"child": 0, "parent": 0}

        if not args.random_off:
            random.shuffle(vcf_files)

        for vcf_file in vcf_files:
            with open(vcf_file, 'r') as vcf_input_file:
                vcf_reader = vcf.Reader(vcf_input_file)
                sampleid = vcf_reader.samples[0]

                if 'runid' in vcf_reader.metadata:
                    runid = vcf_reader.metadata['runid'][0]
                else:  # get runid from sample name (legacy version)
                    runid = "_".join(vcf_file.split("/")[-1].split("bam")[1].split("_")[1:-2])

                if sampleid in merge_dic and runid in merge_dic[sampleid]:
                    """ Check if sample in specific run is merge sample. If so, exclude """
                    excluded_samples_file.write("Sample {sampleid} run {runid} is excluded being merge sample\n".format(
                        sampleid=sampleid,
                        runid=runid)
                    )
                    continue

                if "giab" in sampleid.lower() or "control" in sampleid.lower():
                    """ Remove GIAB and Control samples as these should not be included in the results """
                    excluded_samples_file.write((
                        "Sample {sampleid} run {runid} is excluded being GIAB or Control sample\n"
                    ).format(
                        sampleid=sampleid,
                        runid=runid)
                    )
                    continue

                if sampleid not in vcf_file:
                    excluded_samples_file.write((
                        "Sample {sampleid} form run {runid} does not have the same ID "
                        "in VCF file {vcf_file} and is excluded\n"
                        ).format(
                            sampleid=sampleid,
                            runid=runid,
                            vcf_file=vcf_file
                        )
                    )
                    continue

                if len(vcf_reader.samples) > 1:  # multisample VCF
                    excluded_samples_file.write(
                        "VCF {vcf_file} from run {runid} is excluded because vcf was a multisample VCF\n".format(
                            vcf_file=vcf_file,
                            runid=runid
                        )
                    )
                    continue

                """ Determine Child or Parent status based on sampleid. """
                sampletype = get_familystatus_sampleid(sampleid)
                if sampletype == "unknown":
                    sampletype = get_familystatus_clarity(sampleid, lims)
                if sampletype == "unknown":
                    excluded_samples_file.write(
                        "Sample {sampleid} run {runid} is excluded as no child or parent status could be determined\n".format(
                            sampleid=sampleid,
                            runid=runid
                        )
                    )
                    continue

                count_familystatus[sampletype] += 1

                if 'gender_refset' in vcf_reader.metadata:
                    gender = vcf_reader.metadata['gender_refset'][0]
                    refset = vcf_reader.metadata['exomedepth_reference'][0]
                else:  # get gender and refset based on legacy
                    edreference = vcf_reader.metadata['EDreference'][0]
                    gender = edreference.split("_")[1]
                    refset = edreference.split("_")[2]

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
                    correl = float(record.genotype(sampleid)['CR'])
                    deldupratio = float(record.genotype(sampleid)['PD'])
                    totalcalls = int(record.genotype(sampleid)['TC'])
                    refsamples = int(record.genotype(sampleid)['RS'])

                    if "SM" in record.FORMAT:  # New Franklin format
                        ratio = float(record.genotype(sampleid)['SM'])
                        ccn = float(record.genotype(sampleid)['CN'])
                    else:  # Old Alissa format
                        ratio = float(record.genotype(sampleid)['RT'])
                        ccn = float(record.genotype(sampleid)['CCN'])

                    bf = float(record.genotype(sampleid)['BF'])

                    if(
                        totalcalls < args.totalcallsqc_min or
                        totalcalls > args.totalcallsqc_max or
                        correl < args.correlqc or
                        deldupratio < args.deldupratioqc_min or
                        deldupratio > args.deldupratioqc_max
                    ):
                        excluded_samples_file.write(
                            "Sample {sampleid} run {runid} is excluded because not all QC are above threshold\n".format(
                                sampleid=sampleid,
                                runid=runid
                            )
                        )
                        count_familystatus[sampletype] -= 1
                        break

                    all_event_file.write((
                        "{sampleid}\t{chrom}\t{start}\t{stop}\t"
                        "{gender}\t{refset}\t{calltype}\t{ntargets}\t"
                        "{svlen}\t{ratio}\t{ccn}\t{bf}\t{correl}\t"
                        "{deldupratio}\t{totalcalls}\t{refsamples}\n").format(
                            sampleid=sampleid, chrom=chrom, start=start, stop=stop,
                            gender=gender, refset=refset, calltype=calltype, ntargets=ntargets,
                            svlen=svlen, ratio=ratio, ccn=ccn, bf=bf, correl=correl,
                            deldupratio=deldupratio, totalcalls=totalcalls, refsamples=refsamples
                    ))

                    event = "{chrom}_{start}_{stop}_{calltype}_{ntargets}".format(
                        chrom=chrom,
                        start=start,
                        stop=stop,
                        calltype=calltype,
                        ntargets=ntargets
                    )

                    if event not in event_dic:
                        event_dic[event] = {
                            "parent": {
                                "count": 0,
                                "bf": [],
                                "ratio": [],
                                "correlation": [],
                                "deldupratio": [],
                                "totalcalls": [],
                                "gender": [],
                                "refset": []
                                },
                            "child": {
                                "count": 0,
                                "bf": [],
                                "ratio": [],
                                "correlation": [],
                                "deldupratio": [],
                                "totalcalls": [],
                                "gender": [],
                                "refset": []
                                }
                            }

                    event_dic[event][sampletype]["count"] += 1
                    event_dic[event][sampletype]["bf"].append(bf)
                    event_dic[event][sampletype]["ratio"].append(ratio)
                    event_dic[event][sampletype]["correlation"].append(correl)
                    event_dic[event][sampletype]["deldupratio"].append(deldupratio)
                    event_dic[event][sampletype]["totalcalls"].append(totalcalls)
                    event_dic[event][sampletype]["gender"].append(gender)  # Not used in output at the moment.
                    event_dic[event][sampletype]["refset"].append(refset)
    return event_dic, count_familystatus["child"], count_familystatus["parent"]


def calculate_statistics(data, status):
    totalcount = data['count']
    stats = {
        'total_{0}'.format(status): totalcount,
        'correlation_avr_{0}'.format(status): "n/a",
        'correlation_stdev_{0}'.format(status): "n/a",
        'deldup_avr_{0}'.format(status): "n/a",
        'deldup_stdev_{0}'.format(status): "n/a",
        'totalcalls_avr_{0}'.format(status): "n/a",
        'totalcalls_stdev_{0}'.format(status): "n/a",
        'bf_avr_{0}'.format(status): "n/a",
        'bf_stdev_{0}'.format(status): "n/a",
        'ratio_avr_{0}'.format(status): "n/a",
        'ratio_stdev_{0}'.format(status): "n/a",
        'bfs_{0}'.format(status): "n/a",
        'ratios_{0}'.format(status): "n/a",
        'refsets_{0}'.format(status): "n/a"
        }

    if totalcount >= 1:  # 1 or more element(s) needed to calculate mean
        """ Subselection of first 20 BF, Ratio, Refset. Sort on BF and sort Ratio and Refset accordingly. """
        bfs_subset, ratios_subset, refset_subset = zip(
            *sorted(zip(data["bf"][0:20], data["ratio"][0:20], data["refset"][0:20]))
        )

        bfs = "<br>{}</br>".format('</br><br>'.join(['%.1f' % elem for elem in bfs_subset]))
        ratios = "<br>{}</br>".format('</br><br>'.join(['%.1f' % elem for elem in ratios_subset]))
        refsets = "<br>{}</br>".format('</br><br>'.join(refset_subset))

        stats["bfs_{0}".format(status)] = bfs
        stats["ratios_{0}".format(status)] = ratios
        stats["refsets_{0}".format(status)] = refsets

        stats["correlation_avr_{0}".format(status)] = "%.3f" % (statistics.mean(data['correlation']))
        stats["deldup_avr_{0}".format(status)] = "%.0f" % (statistics.mean(data['deldupratio']))
        stats["totalcalls_avr_{0}".format(status)] = "%.0f" % (statistics.mean(data['totalcalls']))
        stats["bf_avr_{0}".format(status)] = "%.1f" % (statistics.mean(data['bf']))
        stats["ratio_avr_{0}".format(status)] = "%.2f" % (statistics.mean(data['ratio']))

        if totalcount > 1:  # 2 or more elements neede to calculate stdev
            stats["corellation_stdev_{0}".format(status)] = "%.3f" % (statistics.stdev(data['correlation']))
            stats["deldup_stdev_{0}".format(status)] = "%.0f" % (statistics.stdev(data['deldupratio']))
            stats["totalcalls_stdev_{0}".format(status)] = "%.0f" % (statistics.stdev(data['totalcalls']))
            stats["bf_stdev_{0}".format(status)] = "%.1f" % (statistics.stdev(data['bf']))
            stats["ratio_stdev_{0}".format(status)] = "%.2f" % (statistics.stdev(data['ratio']))
    return stats


def make_bed_detail(args, event_dic, children, parents):
    event_file_ucsc = "{outputfolder}/{outputfile}_UCSC.bed".format(
        outputfolder=args.outputfolder, outputfile=args.outputfile
    )
    event_file_igv = "{outputfolder}/{outputfile}_IGV.bed".format(
        outputfolder=args.outputfolder, outputfile=args.outputfile
    )
    event_file_alissa = "{outputfolder}/{outputfile}_Alissa.bed".format(
        outputfolder=args.outputfolder, outputfile=args.outputfile
    )
    with \
            open(event_file_ucsc, "w") as event_file_ucsc, \
            open(event_file_igv, "w") as event_file_igv, \
            open(event_file_alissa, "w") as event_file_alissa:

        """ Write header in BED file """
        event_file_ucsc.write((
            "track name=\"{outputfile}\" type=\"bedDetail\" description=\"CNVs called by Exomedepth "
            "using {outputfile} callset. #Child={children} #Parents={parents} \" visibility=3 itemRgb=\"On\"\n"
            ).format(
                outputfile=args.outputfile,
                children=children,
                parents=parents
            )
        )
        event_file_igv.write((
            "track name=\"{outputfile}\" type=\"bed\" description=\"CNVs called by Exomedepth using {outputfile} callset. "
            "#Child={children} #Parents={parents} \" visibility=3 itemRgb=\"On\"\n"
            ).format(
                outputfile=args.outputfile,
                children=children,
                parents=parents
            )
        )
        total_events = []
        alissa_events = []
        for cnv_event in event_dic:
            chrom, start, stop, calltype, ntargets = cnv_event.split("_")
            """ Make start 0-based for BED file """
            alissa_event = [chrom, start, stop]  # needs to be 1-based
            start = int(start) - 1
            event = [chrom, start, stop]  # Add fields Chromosome, start, stop
            total_count = int(event_dic[cnv_event]["parent"]["count"]) + int(event_dic[cnv_event]["child"]["count"])
            all_counts_ratio = event_dic[cnv_event]["parent"]["ratio"] + event_dic[cnv_event]["child"]["ratio"]
            all_counts_bf = event_dic[cnv_event]["parent"]["bf"] + event_dic[cnv_event]["child"]["bf"]
            median_ratio = statistics.median(all_counts_ratio)
            median_bf = statistics.median(all_counts_bf)

            summary = "{total_count}x:RT={median_ratio:.2f}:BF={median_bf:.0f}:NT={ntargets}".format(
                total_count=total_count,
                median_ratio=median_ratio,
                median_bf=median_bf,
                ntargets=ntargets
            )

            event.append(summary)  # Add field 'name'

            """ Append unused but necessary columns for bed file """
            event.extend([0, '.', start, stop])  # Add fields Score, Strand, thickStart, thickEnd

            """ Append color code based on DEL (red) or DUP (blue) """
            if median_ratio >= 1:
                color = "0,0,255"  # Add field itemRgb for duplication
                alissa_event.append('DUP')
            else:
                color = "255,0,0"  # Add field itemRgb for deletion
                alissa_event.append('DEL')
            event.append(color)

            # Append blockCount, blockSizes, blockStarts.
            # Note that these are not used for ntargets as start/stop for exons is unknown.
            event.extend([1, int(stop) - int(start), 0])  # Add fields blockCount, blockSizes, and blockStarts

            #  Append custom annotation field as html in BEDdetail format
            #  Calculate statistics for children.
            substitute_dic = calculate_statistics(event_dic[cnv_event]["child"], 'child')
            #  Calculate statistics for parents.
            substitute_dic.update(calculate_statistics(event_dic[cnv_event]["parent"], 'parent'))

            template_file = Template(open(settings.html).read())
            new_file = template_file.substitute(substitute_dic)
            event.append("".join(new_file.split("\n")))
            total_events.append(event)
            alissa_event.append(total_count)
            alissa_events.append(alissa_event)

        total_events.sort(key=lambda pos: (settings.chromosome_order[pos[0]], int(pos[1]), int(pos[2])))
        alissa_events.sort(key=lambda pos: (settings.chromosome_order[pos[0]], int(pos[1]), int(pos[2])))

        for event in total_events:
            event[0] = "chr{original}".format(original=event[0])
            joined_event = "\t".join([str(i) for i in event])
            event_file_ucsc.write("{joined}\n".format(joined=joined_event))

            joined_event_igv = "\t".join([str(i) for i in event[0:-1]])
            event_file_igv.write("{joined}\n".format(joined=joined_event_igv))

        for alissa_event in alissa_events:
            joined_alissa_event = "\t".join([str(i) for i in alissa_event])
            event_file_alissa.write("{joined}\n".format(joined=joined_alissa_event))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfolder', help='Path to folder including VCF files')
    parser.add_argument('outputfolder', help='Path to output folder')
    parser.add_argument('outputfile', help='output prefix filename')
    parser.add_argument(
        'merge_samples',
        help='Path to file including all merge samples (tab delimited file with column 1 = sampleID, 2 = run/projectID)'
    )
    parser.add_argument(
        '--totalcallsqc_min', default=settings.number_calls[0], type=int,
        help='Threshold for minimum allowed CNV calls (default = 35)'
    )
    parser.add_argument(
        '--totalcallsqc_max', default=settings.number_calls[1], type=int,
        help='Threshold for maximum allowed CNV calls (default = 200)'
    )
    parser.add_argument(
        '--deldupratioqc_min', default=settings.del_dup_ratio[0], type=float,
        help='Threshold for minimum allowed deldupratio (default = 15)'
    )
    parser.add_argument(
        '--deldupratioqc_max', default=settings.del_dup_ratio[1], type=float,
        help='Threshold for maximum allowed deldupratio (default = 85)'
    )
    parser.add_argument(
        '--correlqc', default=settings.correlation, type=float,
        help='Threshold for minimum allowed correlation score (default = 0.98)'
    )
    parser.add_argument(
        '--random_off', action='store_true',
        help=(
            'Do not apply randomisation in VCF-file order.'
            ' This is usefull for validation purposes if same results are needed (default = random selection on)'
        )
    )
    args = parser.parse_args()

    lims = genologics.lims.Lims(settings.clarity_baseuri, settings.clarity_username, settings.clarity_password)
    genologics.lims.TIMEOUT = settings.api_timeout

    if not os.path.isdir(args.outputfolder):
        os.makedirs("{outputfolder}".format(outputfolder=args.outputfolder), exist_ok=True)
    merge_dic = make_merge_dic(args.merge_samples)
    event_dic, children, parents = slice_vcf(args, merge_dic, lims)
    make_bed_detail(args, event_dic, children, parents)
