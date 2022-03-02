#! /usr/bin/env python
import sys

bed_file = open(sys.argv[1], "r")  # BED file
tlen = float(sys.argv[2])  # targeted length of target/exons

for line in bed_file:
    splitline = line.split()
    leng = float(splitline[2]) - float(splitline[1])
    if leng >= (2 * tlen):
        regions = round(leng / tlen)
        chunk = round(leng / regions)
        start = int(splitline[1])
        for i in range(0, regions - 1):
            stop = start + chunk
            print("{chr}\t{start}\t{stop}".format(chr=splitline[0], start=start, stop=stop))
            start += chunk
        print("{chr}\t{start}\t{stop}".format(chr=splitline[0], start=start, stop=splitline[2]))
    else:
        print(line.rstrip())

bed_file.close()
