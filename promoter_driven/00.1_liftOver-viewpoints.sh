WORKDIR=/well/mccarthy/users/jason/projects/wtsa/promoter_driven
INDIR=$WORKDIR/input_files
INFILE=$INDIR/promoter_centric_viewpoints_hg19_adj-RUNX-coords.txt
# Note the input file must have "chr" strings before chromosome numbers
# Also note that I needed to manually swap the start and end coordinates for RUNX_A and RUNX_B
# as liftover threw out an error regarding unacceptable orientation (i.e. end is less than start)
LIFTOVER=/well/mccarthy/users/jason/software/liftOver/liftOver
CHAINFILE=/well/mccarthy/users/jason/software/liftOver/chain_files/hg19ToHg38.over.chain.gz

# Generate input files
cat $INFILE | cut -f 2,3,4 > $INDIR/temp1.bed
cat $INFILE | cut -f 5,6,7 > $INDIR/temp2.bed
cat $INFILE | cut -f 1 > $INDIR/col1.txt
cat $INFILE | cut -f 8 > $INDIR/col8.txt
cat $INFILE | cut -f 9 > $INDIR/col9.txt

$LIFTOVER $INDIR/temp1.bed $CHAINFILE $INDIR/temp1.38.bed $INDIR/unmapped_file1
$LIFTOVER $INDIR/temp2.bed $CHAINFILE $INDIR/temp2.38.bed $INDIR/unmapped_file2

paste $INDIR/col1.txt $INDIR/temp1.38.bed $INDIR/temp2.38.bed $INDIR/col8.txt $INDIR/col9.txt > $INDIR/promoter_centric_viewpoints_hg38.txt

rm $INDIR/temp1.38.bed
rm $INDIR/temp2.38.bed
rm $INDIR/col1.txt
rm $INDIR/col8.txt
rm $INDIR/col9.txt
rm $INDIR/unmapped_file1
rm $INDIR/unmapped_file2
rm $INDIR/temp1.bed
rm $INDIR/temp2.bed

# Note: you must remove "chr" strings from the output file 
