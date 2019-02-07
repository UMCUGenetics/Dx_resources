#!/bin/bash

# Simple script to determine average DP per SNP in current folder

design_vcf='/hpc/cog_bioinf/diagnostiek/production/Dx_tracks/fingerprint/81SNP_design.vcf'
SnpSift='java -Xmx1G -jar /hpc/local/CentOS7/cog_bioinf/snpEff_v4.2/SnpSift.jar'
snp_positions=`cat $design_vcf | grep -v '#' | cut -f2`
samples=`ls -1 *.vcf | wc -l`

for snp in $snp_positions; do
	if [[ $snp != '2847439' && $snp != '2847910' ]]
		then
			printf "$snp\t"
            dp_sum=`cat *.vcf | $SnpSift filter POS=$snp | $SnpSift extractFields - DP | awk '{sum+=$1} END {print sum}'`
            dp_avg=`expr $dp_sum / $samples`
			echo $dp_avg
		else
			printf "$SNP\t"
            dp_avg=`cat *.vcf | $SnpSift filter POS=$snp | $SnpSift filter 'isRef(GEN[0])' | $SnpSift extractFields - DP | grep -v "#" | awk '{sum+=$1; if ($1) n++; } END {print sum / n}'`
			echo $dp_avg
	fi
done
