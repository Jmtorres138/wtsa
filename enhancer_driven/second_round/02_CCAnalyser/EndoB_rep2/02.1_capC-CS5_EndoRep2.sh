#!/bin/bash

###########################################
# The pipeline is located here :
pipePath="/well/mccarthy/users/jason/projects/wtsa/software/CCseqBasicS-VS1.0.8/"
###########################################
pwd

# The fastq files are located here (these cannot be .gz packed files)

Read1="/well/mccarthy/users/jason/projects/wtsa/enhancer_driven/second_round/fastq_files/original_sequencing/Endo_B_R1.fastq.gz"
Read2="/well/mccarthy/users/jason/projects/wtsa/enhancer_driven/second_round/fastq_files/original_sequencing/Endo_B_R2.fastq.gz"

# This is the genome build for the run :

Genome="hg38"


# This is the oligo coordinate input file (containing the dpnII fragments within which your biotinylated capture oligos are) :
# Make this file with this tutorial : http://sara.molbiol.ox.ac.uk/public/telenius/captureManual/Generating_Oligo_Coordinate_File_For_CaptureC_analyser.pdf

OligoFile="/well/mccarthy/users/jason/projects/wtsa/enhancer_driven/second_round/input_files/Oligo-file_merged_GRCh38_no-chr.txt"


# This is the path to the public folder we want to put the data into (this folder can exist, but it does not need to) :

PublicPath="/well/mccarthy/users/jason/projects/wtsa/enhancer_driven/second_round/public/"

# This is the name we give our sample (no fancy characters, don't start the name with a number)

Sample="EndoB_rep2"


Blatpath="/well/mccarthy/users/jason/projects/wtsa/enhancer_driven/second_round/01_REUSE_blat/"

#############################################
# We are now here (the current directory) - the directory we are now running the pipe :

rundir=$( pwd )


# Tell where we are - print it to output :

echo "Running pipe.sh in folder ${rundir}"


# Tell the command to user - print it to output :
echo "${pipePath}/CCseqBasic5.sh.sh -o "${OligoFile}" -s "${Sample}" --pf "${PublicPath}" --genome ${Genome} --chunkmb 1012 --R1 ${Read1} --R2 ${Read2}  --gz --BLATforREUSEfolderPath ${Blatpath} --onlyCCanalyser"


# Run the command :

${pipePath}/CCseqBasic5.sh -o "${OligoFile}" -s "${Sample}" --pf "${PublicPath}" --genome ${Genome} --chunkmb 1012 --R1 ${Read1} --R2 ${Read2} --gz  --BLATforREUSEfolderPath ${Blatpath} --onlyCCanalyser



# Tell the user that we are now finished :

echo "All done !"
echo
