#!/bin/bash

# Simple script to determine average coverage per SNP in current folder
# Would like to add chromosome numbers (column 1 in 81_snps_mip_design_nijmegen_sort.vcf) to output

a=`cat /hpc/cog_bioinf/diagnostiek/production/Dx_resources/fingerprint/81_snps_mip_design_nijmegen_sort.vcf | grep -v '#' | cut -f2`

for SNP in $a; do
	if [[ $SNP != '2847439' && $SNP != '2847910' ]] 
		then
			printf "$SNP\t"
			b=`cat *vcf | grep $SNP | grep 0/0 | cut -f10 | sed 's/:/\t/g' | cut -f2 | awk '{sum+=$1} END {print sum}'`
			c=`cat *vcf | grep $SNP | egrep '0/1|1/1' | cut -f10 | sed 's/:/\t/g' | cut -f3 | awk '{sum+=$1} END {print sum}'`
			d=`cat *vcf | grep $SNP | wc -l`
			e=`expr $b + $c`
			f=`expr $e / $d`
			echo $f
		else
			printf "$SNP\t"
			g=`cat *vcf | grep $SNP | grep 0/0 | cut -f10 | sed 's/:/\t/g' | cut -f2 | awk '{sum+=$1} END {print sum}'`
			h=`expr $g / $d`
			echo $h
	fi
done	
