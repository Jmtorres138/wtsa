#!/bin/bash
#$ -cwd
#$ -P mccarthy.prjc
#$ -N endoRep1_ccanalyser
#$ -q long.qc
#$ -o /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/second_round/logs/endoRep1_ccanalyser.out
#$ -e /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/second_round/logs/endoRep1_ccanalyser.err
#$ -pe shmem 16

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

sh /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/second_round/02_CCAnalyser/EndoB_rep1/02.1_capC-CS5_EndoRep1.sh
