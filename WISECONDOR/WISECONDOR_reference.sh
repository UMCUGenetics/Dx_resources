module load python/2.7.10

## Excecute this script in the folder with all (softlinked) BAMs. 
## BAMs should be located within one folder.
## Disabled: Argument $1 is outputname for reference file, $2 outputname for GCcount reference 

#REF_OUT=$1
#GCC_OUT=$2
GCC_OUT="GCCcount_IAP_NIPT"
REF_OUT="Reference_IAP_NIPT"

WC_path="/hpc/local/CentOS7/cog_bioinf/WISECONDOR/bin/"
REF_GENOME="/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta"
REF_DIR="$PWD/dataFiles/"
GIT_DIR="/hpc/cog_bioinf/data/mapping/diagnostiek/DEV_Dx_resources/WISECONDOR/"

mkdir $REF_DIR

echo "Making Reference dataset for WISECONDOR"
echo "WISECONDOR path = "$WC_path
echo "Reference outputfile = " $REF_DIR/$REF_OUT
echo "GCC outputfile= "$REF_DIR/$GCC_OUT


# Make new reference genome GCcount
if [ -s $REF_DIR/$GCC_OUT ]
then
	echo $REF_OUT" exists, skipping"
else
	echo "python $WC_path/countgc.py $REF_GENOME $REF_DIR/$GCC_OUT"
        python $WC_path/countgc.py $REF_GENOME $REF_DIR/$GCC_OUT
fi

# Pickle bam files
echo "$GIT_DIR/Pickle_bams.py $WC_path $PWD/"
$GIT_DIR/Pickle_bams.py $WC_path $PWD/ 

# GCC bam files
echo "$GIT_DIR/Gcc_bams.py $WC_path $PWD/ $REF_DIR/$GCC" 
$GIT_DIR/Gcc_bams.py $WC_path $PWD/ $REF_DIR/$GCC_OUT

#Make new reference
echo "python $WC_path/newref.py $PWD $REF_DIR/$REF_OUT"
python $WC_path/newref.py $PWD $REF_DIR/$REF_OUT
