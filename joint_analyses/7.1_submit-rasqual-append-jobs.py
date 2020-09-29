#!/usr/bin/python -O
# Jason Matthew Torres
'''
module load Python/3.7.4-GCCcore-8.3.0
python 07.1_submit-rasqual-append-jobs.py
'''
# libraries
import sys,os
import subprocess as sp

work_dir = "/well/mccarthy/users/jason/projects/wtsa/joint_analyses/"
coloc_dir = work_dir + "coloc_analysis_files/"
rasqual_dir = coloc_dir + "rasqual_output/"
job_dir = coloc_dir + "jobs/caqtl/"
log_dir = coloc_dir +"logs/caqtl/"

def submit_append_job(chromo):
    script = '''#!/bin/bash
#$ -P mccarthy.prjc
#$ -N append-job-%s
#$ -q long.qc
#$ -pe shmem 2
#$ -o %s_append.out
#$ -e %s_append.err

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

module load R/3.6.2-foss-2019b
Rscript %s %s
    ''' % (chromo,log_dir+"chr"+chromo,log_dir+"chr"+chromo,
    work_dir+"7.0_rasqual-caqtls_add-header-pvals.R",chromo)
    fout = open(job_dir + "chr" + chromo + "_append.job.sh",'w')
    fout.write(script)
    fout.close()
    command = ["qsub",job_dir+"chr" + chromo+"_append.job.sh"]
    #command = ["sh",job_dir+"chr" + chromo+"_append.job.sh"] # run serially
    sp.check_call(command)

def main():
    chrom_list = ["2","3","4","5","7","8",
                  "9","10","11"]
    #chrom_list = ["18"]
    for chrom in chrom_list:
        fname = rasqual_dir + "with_p_values/results_rasqual_chr" +  chrom + "_with_Pvalues.txt.gz"
        if os.path.exists(fname)==False:
            submit_append_job(chrom)
            #print(fname)

if (__name__=="__main__"):
     main()
