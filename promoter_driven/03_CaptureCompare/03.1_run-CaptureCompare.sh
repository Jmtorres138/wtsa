#!/bin/bash
#$ -cwd

############################
###     Input Options    ###
############################

# Choice of comparison (either "2way" or "3way")
analysis="2way"

# Directory containing 6/9 run files with structure described in README.md.
##path="/well/mccarthy/users/jason/projects/wtsa/promoter_driven/02_CCAnalyser/"
path="/well/mccarthy/users/jason/projects/wtsa/promoter_driven/02_CCAnalyser/"


# Desired output name of samples with the exact name of the respective input folders.
sample1="EndoB"
sample1Directories="EndoB_rep1,EndoB_rep2,EndoB_rep3"

sample2="LCL"
sample2Directories="LCL_rep1,LCL_rep2,LCL_rep3"

sample3=""                                             #if running 2way leave "" blank
sample3Directories=""             #if running 2way leave "" blank

# Name for the analysis run: e.g. "GWAS_Promoter"
name="T2D_Promoter"

# Genome (supports "hg18","hg19","mm9","mm10")
genome="hg38"

# CaptureC analyser version (Supports any version with output strcuture: F6_greenGraphs_combined_Samples_Version format) e.g. "CS5"
version="CS5"

# Path to file containing viewpoints and windowing parameters -
# Format: Viewpoint    Chr Frag_start  Frag_stop Exlusion_Start Exlusion_Stop Plot_Region_Start Plot_Region_Stop Bin_Size Window_size
parameters="/well/mccarthy/users/jason/projects/wtsa/promoter_driven/03_CaptureCompare/capture-compare-input_hg38.txt"

# Path to where you would like the public UCSC hub to be made
public="/well/mccarthy/users/jason/projects/wtsa/promoter_driven/public "

# Name of enzyme used: supports, "dpnII", "nlaIII", "hindIII"
enzyme="dpnII"

#====================================================#
### Core paths for files and scripts.
#====================================================#

#Path to cis3way shell.
CompareShell="/well/mccarthy/users/jason/projects/wtsa/software/CaptureCompare-master/captureCompare/cis3way.sh"

#Default annotation for plots is ref seq genes but any bed file can be inserted.
annotation="/well/mccarthy/users/jason/projects/wtsa/promoter_driven/input_files/reference_genes/RefSeqGenes_${genome}.bed"

#Path to folder containing digested genomes - will search here for correct digest and create if not present.
digest="/well/mccarthy/users/jason/projects/wtsa/digests/"

#====================================================#
### Run Analysis
#====================================================#

# We are now here
rundir=$( pwd )
echo "Running ${CompareShell} in ${rundir}"

# Print run command
echo "${CompareShell} --analysis=${analysis} --path=${path} --sample1=${sample1} --sample2=${sample2} --sample3=${sample3} --directories=${sample1Directories},${sample2Directories},${sample3Directories} --name=${name} --genome=${genome} --version=${version} --parameters=${parameters} --annotation=${annotation} --public=${public} --restrictionenzyme=${enzyme} --frags=${digest}"

# Run the command :
${CompareShell} --analysis=${analysis} --path=${path} --sample1=${sample1} --sample2=${sample2} --sample3=${sample3} --directories=${sample1Directories},${sample2Directories},${sample3Directories} --name=${name} --genome=${genome} --version=${version} --parameters=${parameters} --annotation=${annotation} --public=${public} --restrictionenzyme=${enzyme} --frags=${digest}

echo "All done !"
echo
