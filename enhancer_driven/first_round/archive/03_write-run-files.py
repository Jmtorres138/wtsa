import sys,os

cur_dir = "/t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/first_round/"

oligofile_path = cur_dir + "Oligofile.txt"
pub_path = "/public/jtorres/wtsa/enhancer_driven/first_round/"
blatdir_path = cur_dir + "REUSE_blat_renamed"
samp_list = ["Endo_A","Endo_B","Endo_C","hESC_A","hESC_B","hESC_C"]


for samp in samp_list:
    samp_dir = cur_dir + samp + "/"
    if os.path.exists(samp_dir) == True:
        fout = open(samp_dir + "run.sh",'w')
        script = '''#!/bin/bash

###########################################
# The pipeline is located here :
pipePath="/home/molhaem2/telenius/CCseqBasic/CM5"
###########################################
# This is the name we give our sample (no fancy characters, don't start the name with a number)

Sample="%s"

#############################################

# The fastq files are located here (these cannot be .gz packed files)

# This is the genome build for the run :

Genome="hg19"

# This is the oligo coordinate input file (containing the dpnII fragments within which your biotinylated capture oligos are) :
# Make this file with this tutorial : http://sara.molbiol.ox.ac.uk/public/telenius/captureManual/Generating_Oligo_Coordinate_File_For_CaptureC_analyser.pdf

OligoFile="%s"

# We are now here (the current directory) - the directory we are now running the pipe :

rundir=$( pwd )

# This is the path to the public folder we want to put the data into (this folder can exist, but it does not need to) :
PublicPath="%s${Sample}"

# Reuse the blat runs from here :

reusedir="%s"

#############################################
# We are now here (the current directory) - the directory we are now running the pipe :

rundir=$( pwd )


# Tell where we are - print it to output :

echo "Running pipe.sh in folder ${rundir}"

# Tell the command to user - print it to output :

echo "${pipePath}/pipeRainbow.sh -o "${OligoFile}" -s "${Sample}" --pf "${PublicPath}" --genome ${Genome} --chunkmb 1012 --BLATforREUSEfolderPath ${reusedir} --useSymbolicLinks --useClusterDiskArea --wholenode24 --bowtie1 -v 3 --wobblyEndBinWidth 20"


# Run the command :


${pipePath}/pipeRainbow.sh -o "${OligoFile}" -s "${Sample}" --pf "${PublicPath}" --genome ${Genome} --chunkmb 1012 --BLATforREUSEfolderPath ${reusedir} --useSymbolicLinks --useClusterDiskArea --wholenode24 --bowtie1 -v 3 --wobblyEndBinWidth 20


# Tell the user that we are now finished :

echo
date
echo "All done !"
echo
        ''' % (samp,oligofile_path,pub_path,blatdir_path)
        fout.write(script+"\n")
        fout.close()
    else:
        print("Directory does not exist: Run and check output of 02 script")
