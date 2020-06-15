#! /usr/bin/env python
import sys

bed_file = open(sys.argv[1],"r")  # BED file
tlen = int(sys.argv[2]) # targeted length of target/exons

for line in bed_file:
    splitline = line.split()
    leng = int(splitline[2]) - int(splitline[1])
    if leng >= (2*tlen):
        regions = int(round(float(leng)/float(tlen)))
        chunk = int(round(leng/regions))
        start = int(splitline[1])
        for i in range(0, regions - 1):
           stop = start + chunk
           print splitline[0] + "\t" + str(start) + "\t" + str(stop)
           start += chunk + 1
        print splitline[0] + "\t" + str(start) + "\t" + str(splitline[2])
    else:
        print line.rstrip()
       

