#!/bin/bash
#$ -cwd
#$ -P mccarthy.prjc
#$ -N hescRep1_ccanalyser
#$ -q long.qc
#$ -o /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/second_round/logs/hescRep1_ccanalyser.out
#$ -e /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/second_round/logs/hescRep1_ccanalyser.err
#$ -pe shmem 16

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

sh /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/second_round/02_CCAnalyser/hESC_rep1/02.1_capC-CS5_hescRep1.sh
