#!/bin/bash
folder=$1
bed=$2
mail=$3

for bam in $folder/*realigned.bam ; do
    #echo "/hpc/local/CentOS7/cog_bioinf/sambamba_v0.6.5/sambamba_v0.6.5 depth region $bam -L $bed -m -q 10 -F \"mapping_quality >= 20 and not duplicate and not failed_quality_control and not secondary_alignment\" -t 6 -o $bam\_coverage"|qsub -cwd -pe threaded 6 -l h_vmem=40G -l h_rt=10:00:00 -m a -M m.elferink@umcutrecht.nl 

    sbatch -t 4:0:00 --gres=tmpspace:20G --mem=20G -c 2  --mail-type=FAIL --mail-user=$mail --wrap="/hpc/local/CentOS7/cog_bioinf/sambamba_v0.6.5/sambamba_v0.6.5 depth region $bam -L $bed -m -q 10 -F \"mapping_quality >= 20 and not duplicate and not failed_quality_control and not secondary_alignment\" -t 6 -o $bam\_coverage"

done
