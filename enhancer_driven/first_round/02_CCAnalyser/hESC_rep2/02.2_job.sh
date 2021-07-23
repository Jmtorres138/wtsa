#!/bin/bash
#$ -cwd
#$ -P mccarthy.prjc
#$ -N E1_hescRep2_ccanalyser
#$ -q short.qc
#$ -o /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/first_round/logs/hescRep2_ccanalyser.out
#$ -e /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/first_round/logs/hescRep2_ccanalyser.err
#$ -pe shmem 16

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

sh /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/first_round/02_CCAnalyser/hESC_rep2/01.1_capC-CS5_hescRep2.sh
