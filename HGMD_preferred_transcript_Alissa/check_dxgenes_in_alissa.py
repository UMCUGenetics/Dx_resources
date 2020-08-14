#! /usr/bin/env python
import sys
import os
import re
import glob
import argparse
from datetime import datetime
import datetime
import time

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dxgenes', type=argparse.FileType('r'), help='File containing all relevant diagnostics genes (1 column HGNC geneID)')
    parser.add_argument('alissa', type=argparse.FileType('r'), help='HGMD dump file from Alissa (2 columns: 1) HGNC geneID 2)NM_TranscriptID with version)')
    parser.add_argument('alias', type=argparse.FileType('r'), help='Alias table (id=Alias gene_id= GeneID in Exoncov)')
    args = parser.parse_args()

    dic_alissa = {}
    for line in args.alissa:
        splitline = line.split()
        if len(splitline) == 2:
            if splitline[0] not in dic_alissa:
                dic_alissa[splitline[0]] = [splitline[1]]
            else:
                dic_alissa[splitline[0]] += [splitline[1]]

    dic_alias = {}
    alias_double = {}
    for line in args.alias:
        splitline = line.split()
        if splitline[1] not in dic_alias:
            dic_alias[splitline[1]] = [splitline[0]]
        else:
            dic_alias[splitline[1]] += [splitline[0]]

        if splitline[0] not in alias_double:
            alias_double[splitline[0]] = [splitline[1]]
        else:
            alias_double[splitline[0]] += [splitline[1]]


    print_dic = {"Detected_in_Alissa":[], "Detected_in_Alissa_as_Alias":[], "Not_Detected_in_Alissa":[]}
    for line in args.dxgenes:
        splitline = line.split()
        if splitline[0] in dic_alissa:  #If GeneID is detected in Alissa list, print all gene + transcript
            if len(dic_alissa[splitline[0]]) == 1:
                print_dic["Detected_in_Alissa"] += [[splitline[0], dic_alissa[splitline[0]][0]]]
            else:
                for transcript in dic_alissa[splitline[0]]:
                    print_dic["Detected_in_Alissa"] += [[splitline[0], transcript]]
        else:  #If geneID is not detected, try to recover the alias.
            """Check if Alias is in Alissa"""
            if splitline[0] in dic_alias:  #If gene has alias, print geneid + alias information
                detected = False
                for item in dic_alias[splitline[0]]:
                    if item in dic_alissa:
                        if item in alias_double and len(alias_double[item])>1:
                            if len(dic_alissa[item]) == 1:
                                print_dic["Detected_in_Alissa_as_Alias"] += [[splitline[0],dic_alissa[item][0],item,"WARNING:recursive geneID"]]
                            else:
                                for transcript in dic_alissa[item]:
                                    print_dic["Detected_in_Alissa_as_Alias"] += [[splitline[0],transcript,item,"WARNING:recursive geneID"]]
                        else:
                            if len(dic_alissa[item]) == 1:
                                print_dic["Detected_in_Alissa_as_Alias"] += [[splitline[0],dic_alissa[item][0],item]]
                            else:
                                for transcript in dic_alissa[item]:
                                    print_dic["Detected_in_Alissa_as_Alias"] += [[splitline[0],transcript,item,"WARNING:recursive geneID"]]
                        detected = True
                if detected == False:
                     print_dic["Not_Detected_in_Alissa"] += [[splitline[0]]]
            else:
                print_dic["Not_Detected_in_Alissa"] += [[splitline[0]]]

    print "Remark\tGene(Exoncov)\tTranscript(Alissa)\tAlias(Alissa)"
    for item in print_dic["Detected_in_Alissa"]:
        print "Detected_in_Alissa\t","\t".join(item)
    for item in print_dic["Detected_in_Alissa_as_Alias"]:
        print "Detected_in_Alissa_as_Alias\t","\t".join(item)
    for item in print_dic["Not_Detected_in_Alissa"]:
        print "Not_Detected_in_Alissa\t","\t".join(item)
