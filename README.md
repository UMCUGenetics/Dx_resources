# Repository with standalone wrapper scripts



**get_stats_from_flagstat.pl**:	is used within the cleanup_Dx scripts (EXOME and GENOME) to calculate flagstats statistics.


folder **ExomeDepth**:	contains the wrapper scripts for Exomedepth CNV analysis\
			see ExomeDepth/README.md for more information

folder **ParseCNVQC**: contains script to make summary of CNV QC's and send mail. Note: this tool is not for diagnostic use!

**/GIAB_rtgtools/compare_vcf_to_giab_truth.pl**: is used for GIAB comparisons: see GENOOM160


**mips_fingerprint/coverage_mips.sh**: is used for validation of new MIP pool, see GENOOM080 5.3.3


**HGMD_preferred_transcript_Alissa/check_dxgenes_in_alissa.py**: script to check if Alissa HGMD transcripts are present in the UMCU target BED file (i.e. used in Exoncov)
Usage:\
``` bash
python check_dxgenes_in_alissa.py {gene file} {Alissa HGMD file} {gene alias}  > {output} 
```
{gene file} = list with all genes in the target BED file\
{Alissa HGMD file} = tab seperated file with geneID and transcriptID\
{gene alias} = tab seperated file with Alias_geneID and HGNC_approved_geneID

