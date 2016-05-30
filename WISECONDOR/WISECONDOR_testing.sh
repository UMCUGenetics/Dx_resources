module load python/2.7.10

SAMPLE=$1
WC_path="/hpc/local/CentOS7/cog_bioinf/WISECONDOR/bin/"
REF_gcc=/hpc/cog_bioinf/data/mapping/diagnostiek/DEV_Dx_resources/WISECONDOR/dataFiles/GCCcount_IAP_NIPT
REF_ref=/hpc/cog_bioinf/data/mapping/diagnostiek/DEV_Dx_resources/WISECONDOR/dataFiles/Reference_IAP_NIPT
TEST_folder=$PWD/WISECONDOR/$SAMPLE/
GIT_DIR="/hpc/cog_bioinf/data/mapping/diagnostiek/DEV_Dx_resources/WISECONDOR/"


echo "Sample for testing is "$SAMPLE
echo "WISECONDOR path is "$WC_path

mkdir $PWD/WISECONDOR/
mkdir $PWD/WISECONDOR/$SAMPLE
cd $PWD/WISECONDOR/$SAMPLE
ln -sd $PWD/../../$SAMPLE/mapping/*dedup.ba* .

echo "$GIT_DIR/Pickle_bams.py $WC_path $TEST_folder"
$GIT_DIR/Pickle_bams.py $WC_path $TEST_folder

echo "$GIT_DIR/Gcc_bams.py $WC_path $TEST_folder $REF_gcc" 
$GIT_DIR/Gcc_bams.py $WC_path $TEST_folder $REF_gcc

echo "$GIT_DIR/Test_bams.py $WC_path $TEST_folder $REF_ref"
$GIT_DIR/Test_bams.py $WC_path $TEST_folder $REF_ref

echo "$GIT_DIR/Plot_bams.py $WC_path $TEST_folder"
$GIT_DIR/Plot_bams.py $WC_path $TEST_folder 

