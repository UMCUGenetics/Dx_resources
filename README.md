# Within this folder there are files needed for IAP, the coverage tool, and NIPT.
BED files
	* ENSEMBL_UCSC_merged_collapsed_sorted_v2.bed	
		This is the 'master' BED file without flanks. 
		Warning: this file is 0-based
		This file contains all exon from ENSEMBL (v75) and USCS (downloaded  07-08-2014). For more details see /hpc/cog_bioinf/data/Martin/Exomes/ENSEMBL_UCSC_merge_bed/
		This BED file is version 2 in /hpc/cog_bioinf/data/Martin/Exomes/ENSEMBL_UCSC_merge_bed/Version2/
		This file is not used directly at the moment.
	* ENSEMBL_UCSC_merged_collapsed_sorted_v2_20bpflank.bed
		Normal BED file, with 20bp flanks. This file should be used for the Exoncov_v3 tool.
		Warning: this file is 0-based
	* ENSEMBL_UCSC_merged_collapsed_sorted_v2_100bpflank.bed
		Normal BED file, with 100bp flanks. This file can be used as track in IGV/Alamut
		Warning: this file is 0-based
	* ENSEMBL_UCSC_merged_collapsed_sorted_v2_20bpflank.list
		List file with 20bp flanks: first 5 columns of BED with header. This file should be used for mapping stats (TARGET)for picard.
		Warning: this file is 1-based
 	* ENSEMBL_UCSC_merged_collapsed_sorted_v2_100bpflank.list
		List file with 100bp flanks: first 5 columns of BED with header. This file is used for variant calling with the haplotype caller.
		Warning: this file is 1-based

Exoncov_v3 files
	gpanels_10062015.txt
		This file contains all currently used gen-panels for Dx.
		Remark: OWS02(Schisis) is not working due to included noncoding genes/genes not included the BED file 	
	Preferred_transcript_list.txt
		This file contains the manually curated/confirmed link between HGNC, NM and ENST.
		Curation is done by the labspecialists.
		Currently, all preferred transcripts for all gen-enrichments are pooled (duplicates are removed).
		It is possible that a specific NM has different linked ENST-ID. This is the result of different choices of the labspecialists.
		Note: if more than 1 preferred transcript is known within this file, Exoncov_v3 skips these preferred transcripts, and choses the longest NM, 
		or if no NM's are known, the longest ENST.
	NM_ENSEMBL_HGNC.txt 
		This file includes the known link between GENE (HGNC) and all known NM/ENST transcripts
Controle samples:
	control_samples_full_sort.vcf
		This is a VCF with all SNPs + genotypes for the 4 control samples. This is a genome wide assay.


