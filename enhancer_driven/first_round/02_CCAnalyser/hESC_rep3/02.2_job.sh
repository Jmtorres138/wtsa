#!/bin/bash
#$ -cwd
#$ -P mccarthy.prjc
#$ -N E1_hescRep3_ccanalyser
#$ -q long.qc
#$ -o /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/first_round/logs/hescRep3_ccanalyser.out
#$ -e /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/first_round/logs/hescRep3_ccanalyser.err
#$ -pe shmem 30

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

sh /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/first_round/02_CCAnalyser/hESC_rep3/02.1_capC-CS5_hescRep3.sh
