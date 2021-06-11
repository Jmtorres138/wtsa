import sys,os

cur_dir = "/t1-data/user/hugheslab/jtorres/wtsa/promoter_driven/"

oligofile_path = "/t1-data/user/hugheslab/jtorres/wtsa/promoter_driven/oligo-file.txt"
blatdir_path = "/t1-data/user/hugheslab/jtorres/promoter-driven/Processed/Blymph_B/F4_blatPloidyFilteringLog_Blymph_B_CB4/BlatPloidyFilterRun/REUSE_blat"
samp_list = ["BlymphA","BlymphB","BlymphC","EndoBA","EndoBC","EndoBD"]


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
PublicPath="/public/jtorres/wtsa/promoter_driven/${Sample}"

# Reuse the blat runs from here :

reusedir="%s"

#############################################
# We are now here (the current directory) - the directory we are now running the pipe :

rundir=$( pwd )


# Tell where we are - print it to output :

echo "Running pipe.sh in folder ${rundir}"

# Tell the command to user - print it to output :

echo "${pipePath}/pipeRainbow.sh -o "${OligoFile}" -s "${Sample}" --pf "${PublicPath}" --genome ${Genome} --chunkmb 1012 --BLATforREUSEfolderPath ${reusedir} --useSymbolicLinks --useClusterDiskArea --wholenode24 --bowtie2 --wobblyEndBinWidth 20"


# Run the command :

${pipePath}/pipeRainbow.sh -o "${OligoFile}" -s "${Sample}" --pf "${PublicPath}" --genome ${Genome} --chunkmb 1012 --BLATforREUSEfolderPath ${reusedir} --useSymbolicLinks --useClusterDiskArea --wholenode24 --bowtie1 -v 3 --wobblyEndBinWidth 20


# ${pipePath}/pipeRainbow.sh -o "${OligoFile}" -s "${Sample}" --pf "${PublicPath}" --genome ${Genome} --chunkmb 1012 --BLATforREUSEfolderPath ${reusedir} --useSymbolicLinks --useClusterDiskArea --wholenode24 --bowtie2 --wobblyEndBinWidth 20
# Tell the user that we are now finished :

#  ${pipePath}/pipeRainbow.sh -o "${OligoFile}" -s "${Sample}" --pf "${PublicPath}" --genome ${Genome} --chunkmb 1012 --BLATforREUSEfolderPath ${reusedir} --useSymbolicL# inks --bowtie2 -p 8 --onlyHub --wobblyEndBinWidth 20
# Encountered a RAM MEMORY Issue so have to remove --useClusterDiskArea --wholenode24 and add -p 8 flag with --onlyCCanalyser

echo
date
echo "All done !"
echo
        ''' % (samp,oligofile_path,blatdir_path)
        fout.write(script+"\n")
        fout.close()
    else:
        print("Directory does not exist: Run and check output of 02 script")
