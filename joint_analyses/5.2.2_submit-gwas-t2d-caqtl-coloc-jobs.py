#!/usr/bin/python -O
# Jason Matthew Torres
'''
module load Python/3.7.4-GCCcore-8.3.0
python 02.4_plink-glm-gwas.py
'''
# libraries
import sys,os
import subprocess as sp

work_dir = "/well/mccarthy/users/jason/projects/wtsa/joint_analyses/"
coloc_dir = work_dir + "coloc_analysis_files/"
gwas_sig_dir = coloc_dir + "conditioned_gwas_summarystats/"
job_dir = coloc_dir + "jobs/caqtl/"
log_dir = coloc_dir +"logs/caqtl/"


def submit_coloc_job(chromo):
    script = '''#!/bin/bash
#$ -P mccarthy.prjc
#$ -N coloc-job-%s
#$ -q long.qc
#$ -o %s.out
#$ -e %s.err

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

module load R/3.6.2-foss-2019b
Rscript %s %s
    ''' % (chromo,log_dir+chromo,log_dir+chromo,
    work_dir+"5.0_caqtl-coloc-by-chrom.R",chromo)
    fout = open(job_dir + chromo + "_job.sh",'w')
    fout.write(script)
    fout.close()
    command = ["qsub",job_dir+chromo+"_job.sh"]
    sp.check_call(command)

def main():
    file_list = os.listdir(gwas_sig_dir)
    chrom_list = list(set([x.split("-")[-2] for x in file_list]))
    for chrom in chrom_list:
        submit_coloc_job(chrom)

if (__name__=="__main__"):
     main()
