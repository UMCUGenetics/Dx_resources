#! /usr/bin/env python3
import argparse
import statistics

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create ExomeDepth Metric summary file')
    parser.add_argument('exomedepth_logs', type=argparse.FileType('r'), nargs='*', help='Exomedepth log files')
    arguments = parser.parse_args()
    stats_dic = {"CR": [], "PD": [], "TC": []}
    for exomedepth_qc_file in arguments.exomedepth_logs:
        for line in exomedepth_qc_file:
            splitline = line.split()
            correlation = float(splitline[3])
            deldupratio = float(splitline[4])
            totalcount = float(splitline[5])
            warnings = splitline[6:]
            stats_dic["CR"].append(correlation)
            stats_dic["PD"].append(deldupratio)
            stats_dic["TC"].append(totalcount)
            print("{sample};CM={model};REFSET={refset};CR={correl};PD={deldupratio};TC={totalcount}\t{warnings}\r".format(
                 sample=splitline[0],
                 model=splitline[1],
                 refset=splitline[2],
                 correl="%.4f" % correlation,
                 deldupratio="%.2f" % deldupratio,
                 totalcount="%.0f" % totalcount,
                 warnings="\t".join(warnings)
                 )),

    print("\r")
    print("#Average_CR={}\r".format("%.4f" % statistics.mean(stats_dic["CR"]))),
    print("#Average_PD={}\r".format("%.2f" % statistics.mean(stats_dic["PD"]))),
    print("#Average_TC={}\r".format("%.2f" % statistics.mean(stats_dic["TC"]))),
    print("\r")
    print("#Median_CR={}\r".format("%.4f" % statistics.median(stats_dic["CR"]))),
    print("#Median_PD={}\r".format("%.2f" % statistics.median(stats_dic["PD"]))),
    print("#Median_TC={}\r".format("%.2f" % statistics.median(stats_dic["TC"]))),
