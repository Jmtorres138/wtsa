# Concatentate fastq files from original and additional sequencing
WORKDIR=/well/mccarthy/users/jason/projects/wtsa/enhancer_driven/first_round
SEQDIR=$WORKDIR/fastq_files
OUTDIR=$SEQDIR/combined
ORIGDIR=$SEQDIR/original_sequencing
ADDDIR=$SEQDIR/additional_sequencing
mkdir $OUTDIR
echo Endo_A
#zcat $ORIGDIR/Endo_A_R1.fastq.gz $ADDDIR/Endo_A_1.fastq.gz > \
#  $OUTDIR/Endo_A_R1.combined.fastq
gzip $OUTDIR/Endo_A_R1.combined.fastq
#zcat $ORIGDIR/Endo_A_R2.fastq.gz $ADDDIR/Endo_A_2.fastq.gz > \
#  $OUTDIR/Endo_A_R2.combined.fastq
gzip $OUTDIR/Endo_A_R2.combined.fastq
echo Endo_B
#zcat $ORIGDIR/Endo_B_R1.fastq.gz $ADDDIR/Endo_B_1.fastq.gz > \
#  $OUTDIR/Endo_B_R1.combined.fastq
gzip $OUTDIR/Endo_B_R1.combined.fastq
#zcat $ORIGDIR/Endo_B_R2.fastq.gz $ADDDIR/Endo_B_2.fastq.gz > \
#  $OUTDIR/Endo_B_R2.combined.fastq
gzip $OUTDIR/Endo_B_R2.combined.fastq
echo Endo_C
#zcat $ORIGDIR/Endo_C_R1.fastq.gz $ADDDIR/Endo_C_1.fastq.gz > \
#  $OUTDIR/Endo_C_R1.combined.fastq
gzip $OUTDIR/Endo_C_R1.combined.fastq
#zcat $ORIGDIR/Endo_C_R2.fastq.gz $ADDDIR/Endo_C_2.fastq.gz > \
#  $OUTDIR/Endo_C_R2.combined.fastq
gzip $OUTDIR/Endo_C_R2.combined.fastq
echo hESC_A
#zcat $ORIGDIR/hESC_A_R1.fastq.gz $ADDDIR/hESC_A_1.fastq.gz > \
#  $OUTDIR/hESC_A_R1.combined.fastq
gzip $OUTDIR/hESC_A_R1.combined.fastq
#zcat $ORIGDIR/hESC_A_R2.fastq.gz $ADDDIR/hESC_A_2.fastq.gz > \
#  $OUTDIR/hESC_A_R2.combined.fastq
gzip $OUTDIR/hESC_A_R2.combined.fastq
echo hESC_B
zcat $ORIGDIR/hESC_B_R1.fastq.gz $ADDDIR/hESC_B_1.fastq.gz > \
  $OUTDIR/hESC_B_R1.combined.fastq
gzip $OUTDIR/hESC_B_R1.combined.fastq
zcat $ORIGDIR/hESC_B_R2.fastq.gz $ADDDIR/hESC_B_2.fastq.gz > \
  $OUTDIR/hESC_B_R2.combined.fastq
gzip $OUTDIR/hESC_B_R2.combined.fastq
echo hESC_C
zcat $ORIGDIR/hESC_C_R1.fastq.gz $ADDDIR/hESC_C_1.fastq.gz > \
  $OUTDIR/hESC_C_R1.combined.fastq
gzip $OUTDIR/hESC_C_R1.combined.fastq
zcat $ORIGDIR/hESC_C_R2.fastq.gz $ADDDIR/hESC_C_2.fastq.gz > \
  $OUTDIR/hESC_C_R2.combined.fastq
gzip $OUTDIR/hESC_C_R2.combined.fastq
