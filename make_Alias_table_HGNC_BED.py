#! /usr/bin/env python
import sys, os, re
import glob
import commands


# ME 28-09-2018: script for producing Alias table between (recent) HGNC and UCMU Exoncov target files.
# Download HGNC file from https://www.genenames.org/ > Complete dataset download links > Complete HGNC dataset (ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt)
# Downloaded as hgnc_complete_set_[date].txt
# make proteincoding+non-coding RNA file
# > cat hgnc_complete_set_[date].txt| egrep -w '(protein-coding|non-coding)  > hgnc_complete_set_[date]_pc_nc.txt
# Run script
# > ./make_Alias_table_HGNC_BED.py hgnc_complete_set_[date]_pc_nc.txt [TARGET BED file (ENSEMBL_UCSC_merged_collapsed_sorted_v3_20bpflank.bed)] | (read -r; printf "%s\n" "$REPLY"; sort -nk1) > Alias_table_HGNC_BED_[date].txt


def unique(list):
    seen = set()
    seen_add = seen.add
    return [ x for x in list if not (x in seen or seen_add(x))]

def invert_dict(d):
    inverted_dict = {}
    for key in d:
        for item in d[key]:
            if "NM"in item: # only NM transcripts are relevant
                try:
                    inverted_dict[item]
                    inverted_dict[item]+=[key]
                except:
                    inverted_dict[item]=key


    return inverted_dict

hgnc=open(sys.argv[1],"r").readlines()
bed=open(sys.argv[2],"r")


# make dictionary of HGNC file (key=gene, value=transcript(s)
dic={}
for line in hgnc:
    splitline=line.split("\t")
    if "NM" in splitline[23]:   ## only for coding genes
        if "|" in splitline[23]:
            transcript=splitline[23].split("|") 
        else:
            transcript=[splitline[23]]

        for item in transcript:
            stripitem=item.replace("\"","") 		# strip "
            stripitem=stripitem.replace("\"","")	# strip "
            stripitem=stripitem.split(".")[0]		# strip .

            try:
                dic[splitline[1]]
                dic[splitline[1]]+=[stripitem]  
            except:
                dic[splitline[1]]=[stripitem]

# Make dictionairy of BED file (key =gene, values = transcript(s))
bed_dic={}
for line in bed:
    splitline=line.split()
    genes=splitline[4].split(":")
    for gene in genes:
        try:
            bed_dic[gene]
            for item in splitline[6].split(":"):
                bed_dic[gene]+=[item]
            bed_dic[gene]= unique(bed_dic[gene])
        except:
            split=splitline[6].split(":")
            if len(split)==1:
                bed_dic[gene]=[splitline[6]]
            else:
                x=1
                for item in split:
                    if x==1:
                        bed_dic[gene]=[item]
                        x+=1
                    else:
                        bed_dic[gene]+=[item] 
            bed_dic[gene]= unique(bed_dic[gene])

# Check if genes in HGNC file exist in BED file.
# if not, check if transcript matches in BED file, and print corresponding gene.
# otherwise, print gene notpresent (indication that GENE in hgnc is not found in BED)

# invert dictionary to make transcripts key
invert_bed_dic=invert_dict(bed_dic)
invert_dic=invert_dict(dic)

print "GeneID_HGNC\tTranscriptID_BED\tGeneID_BED" 
for item in dic:
    try:
        bed_dic[item]	# if HGNC known, do nothing
    except:
        printlist=[]
        for trans in dic[item]:
            try:
                printlist+=[invert_bed_dic[trans]]
            except:						# if HGNC gene/transcript is not known in BED file (i.e. new gene/transcript)
                printlist+=["notpresent"]
        print item,"\t",",".join(dic[item]),"\t",",".join(unique(printlist))
