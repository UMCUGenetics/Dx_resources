#!/bin/bash
sambamba=$1
folder=$2
bed=$3
mail=$4

for bam in $folder/*.bam ; do
    sbatch -t 4:0:00 --gres=tmpspace:20G --mem=20G -c 2  --mail-type=FAIL --mail-user=$mail --wrap="$sambamba depth region $bam -L $bed -m -q 10 -F \"mapping_quality >= 20 and not duplicate and not failed_quality_control and not secondary_alignment\" -t 6 -o $bam\_coverage"

done
