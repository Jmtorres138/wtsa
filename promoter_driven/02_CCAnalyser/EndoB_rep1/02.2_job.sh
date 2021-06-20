#!/bin/bash
#$ -cwd
#$ -P mccarthy.prjc
#$ -N endoR1_ccanalyser
#$ -q short.qc
#$ -o /well/mccarthy/users/jason/projects/wtsa/promoter_driven/logs/endoR1_ccanalyser.out
#$ -e /well/mccarthy/users/jason/projects/wtsa/promoter_driven/logs/endoR1_ccanalyser.err
#$ -pe shmem 20

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

sh /well/mccarthy/users/jason/projects/wtsa/promoter_driven/02_CCAnalyser/EndoB_rep1/01.1_capC-CS5_EndoRep1.sh
