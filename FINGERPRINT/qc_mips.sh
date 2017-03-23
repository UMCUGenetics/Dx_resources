#!/bin/bash

# Simple script to per sample determine:
# 	Column 2 -> number of MIPs with coverage <15X
#	Column 3 -> number of 1/1 calls
#	Column 4 -> calculation of 0 reads in 1/1 calls as indicator for contamination
#	Column 5 -> male/female according to sample name
#	Column 6 -> number of reads on Y-chromosome

# Variant calling done with UnifiedGenotyper, coverage downsampled to max. 250 reads (default), 0/0 calls show no allele frequency
# Future wish to test no downsampling with -dcov and test adding 1 reads in 0/0 calls to indicator for contamination with HaplotypeCaller

for vcf in $1/*.vcf
do
	a="${vcf##*/}"
	printf "$a\t"

	b=`cat $vcf | grep -v "#" | grep -v "Y" | cut -f10 | fgrep ./. | wc -l`
	c=`cat $vcf | grep -v "#" | grep -v "Y" | cut -f10 | grep 0/0 | sed 's/:/\t/g' | awk '$2<15' | wc -l`
	d=`cat $vcf | grep -v "#" | cut -f10 | egrep '0/1|1/1' | sed 's/:/\t/g' | awk '$3<15' | wc -l`
	e=$((b + c + d))
	printf "$e\t"
	
	f=`cat $vcf | grep 1/1 | wc -l `
	printf "$f\t"
	g=`cat $vcf | grep 1/1 | cut -f10 | sed 's/:/\t/g' | sed 's/,/\t/g' | cut -f2,4 | awk '{total += $1;total2 += $2} END {print total/total2}'`
	printf "$g\t"

	h=`printf "$a" | head -c 9 | tail -c 1`
	printf "$h\t"
	i=`cat $vcf | grep -v '#' | cut -f1,10 | grep Y | sed 's/:/\t/g' | cut -f3 | awk '{SUM += $1} END {print SUM}'`
	echo "$i"
done
