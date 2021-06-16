#!/bin/bash
#$ -cwd
#$ -P mccarthy.prjc
#$ -N LCL_R2_ccanalyser
#$ -q short.qc
#$ -o /well/mccarthy/users/jason/projects/wtsa/promoter_driven/logs/LCL_R2_ccanalyser.out
#$ -e /well/mccarthy/users/jason/projects/wtsa/promoter_driven/logs/LCL_R2_ccanalyser.err
#$ -pe shmem 16

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

sh /well/mccarthy/users/jason/projects/wtsa/promoter_driven/02_CCAnalyser/LCL_rep2/01.1_capC-CS5_LCL_rep2.sh
