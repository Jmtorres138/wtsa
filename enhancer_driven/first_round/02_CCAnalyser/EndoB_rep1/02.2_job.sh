#!/bin/bash
#$ -cwd
#$ -P mccarthy.prjc
#$ -N E1_endoRep1_ccanalyser
#$ -q short.qc
#$ -o /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/first_round/logs/endoRep1_ccanalyser.out
#$ -e /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/first_round/logs/endoRep1_ccanalyser.err
#$ -pe shmem 16

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

sh /well/mccarthy/users/jason/projects/wtsa/enhancer_driven/first_round/02_CCAnalyser/EndoB_rep1/01.1_capC-CS5_EndoRep1.sh
