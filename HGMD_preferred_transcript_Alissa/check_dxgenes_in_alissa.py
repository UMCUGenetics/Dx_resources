#! /usr/bin/env python
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('dxgenes_file', type=argparse.FileType('r'), help='File containing all relevant diagnostics genes (1 column HGNC geneID)')
    parser.add_argument('alissa_file', type=argparse.FileType('r'), help='HGMD dump file from Alissa (2 columns: 1) HGNC geneID 2)NM_TranscriptID with version)')
    parser.add_argument('alias_file', type=argparse.FileType('r'), help='Alias table (id=Alias gene_id= GeneID in Exoncov)')
    args = parser.parse_args()

    alissa_genes = {}
    for line in args.alissa_file:
        splitline = line.split()
        if len(splitline) == 2:
            if splitline[0] not in alissa_genes:
                alissa_genes[splitline[0]] = []
            alissa_genes[splitline[0]].append(splitline[1])

    gene_aliases = {}
    gene_aliases_reversed = {}
    for line in args.alias_file:
        gene_alias, gene = line.strip().split()

        if gene not in gene_aliases:
            gene_aliases[gene] = []
        gene_aliases[gene].append(gene_alias)

        if gene_alias not in gene_aliases_reversed:
            gene_aliases_reversed[gene_alias] = []
        gene_aliases_reversed[gene_alias].append(gene)

    print_dic = {"Detected_in_Alissa":[], "Detected_in_Alissa_as_Alias":[], "Not_Detected_in_Alissa":[]}
    for line in args.dxgenes_file:
        dx_gene = line.split()[0]
        if dx_gene in alissa_genes:  #If GeneID is detected in Alissa list, print all gene + transcript
            for transcript in alissa_genes[dx_gene]:
                print_dic["Detected_in_Alissa"].append([dx_gene, transcript])
        else:  #If geneID is not detected, try to recover the alias.
            """Check if Alias is in Alissa"""
            detected = False
            if dx_gene in gene_aliases:  #If gene has alias, print geneid + alias information
                for alias_gene in gene_aliases[dx_gene]:
                    if alias_gene in alissa_genes:
                        if alias_gene in gene_aliases_reversed and len(gene_aliases_reversed[alias_gene])>1:
                            for transcript in alissa_genes[alias_gene]:
                                print_dic["Detected_in_Alissa_as_Alias"].append([dx_gene, transcript, alias_gene, "WARNING:recursive geneID"])
                        else:
                            for transcript in alissa_genes[alias_gene]:
                                print_dic["Detected_in_Alissa_as_Alias"].append([dx_gene, transcript, alias_gene])
                        detected = True
            if detected == False:
                print_dic["Not_Detected_in_Alissa"].append([dx_gene])

    print "Remark\tGene(Exoncov)\tTranscript(Alissa)\tAlias(Alissa)"
    for item in print_dic["Detected_in_Alissa"]:
        print("Detected_in_Alissa\t{item}".format(item="\t".join(item)))
    for item in print_dic["Detected_in_Alissa_as_Alias"]:
        print("Detected_in_Alissa_as_Alias\t{item}".format(item="\t".join(item)))
    for item in print_dic["Not_Detected_in_Alissa"]:
        print("Not_Detected_in_Alissa\t{item}".format(item="\t".join(item)))
