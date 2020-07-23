#!/usr/bin/python -O
# Jason Matthew Torres
'''
module load Python/3.7.4-GCCcore-8.3.0
python 02.4_plink-glm-gwas.py
'''
# libraries
import sys,os
import subprocess as sp

work_dir = "/well/emberson/users/bjk420/projects/freeze_100k/04_evaluate-gwas-methods/"
job_dir = work_dir+"plink/jobs/"
log_dir = work_dir+"plink/logs/"
in_dir = work_dir+"plink/input_files/"
out_dir = work_dir+"plink/output_files/"
king_dir = work_dir + "plink/king/"
plink = "/apps/well/plink/2.00a-20170724/plink2"

pheno_file = in_dir + "BASE_DIABETES.txt"

pre_list = ["rel-pruned-p354_no-twins-dups",
            "rel-pruned-p177_exclude-first-degree",
            "rel-pruned-p0884_exclude-second-degree",
            "rel-pruned-p0442_exclude-third-degree",
            "rel-pruned-p0221_unrelated"]

cov_list = ["covar-age-sex-pc10","covar-age-sex-pc20",
            "covar-bmi-age-sex-pc10","covar-bmi-age-sex-pc20"]

def submit_plink_glm_job(geno_dir,geno_pre,pheno_file,covar_file,
                         out_dir,out_pre):
    # geno_dir is the directory where the plink input files are locate (.bed,.bim,.fam)
    # geno_pre the file prefix for the plink input file_list
    # the file path to the input phenotype file (will test all phenos in the file)
    # the file path to the covariate file

    script = '''
#!/bin/bash
#$ -cwd
#$ -P emberson.prjc
#$ -N glm-job-%s
#$ -q short.qc
#$ -o %s.out
#$ -e %s.err

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

%s --bfile %s%s \
--glm hide-covar \
--pheno %s --1 \
--covar %s \
--out %s%s
    ''' % (out_pre,log_dir+out_pre,log_dir+out_pre,plink,geno_dir,geno_pre,
    pheno_file,covar_file,out_dir,out_pre)
    fout = open(job_dir + out_pre + "_job.sh",'w')
    fout.write(script)
    fout.close()
    command = ["qsub",job_dir+out_pre+"_job.sh"]
    sp.check_call(command)

def main():
    #geno_dir = "/well/emberson/users/bjk420/projects/freeze_100k/plink_subset/subset3/"
    #geno_pre = 'mcps-freeze100k-gt-hg38_subset3'
    #pheno_file = in_dir + "t2d_status.txt"
    #covar_file = in_dir + "mcps-freeze100k-gt-hg38_subset3.covar-pc10.txt"
    #out_pre = geno_pre + ".no-rel-filt.pc10"
    for pre in pre_list:
        for cov in cov_list:
            geno_dir,geno_pre = king_dir,pre
            covar_file = in_dir+pre+"."+cov+".txt"
            out_pre = pre+"."+cov
            #print(out_pre)
            submit_plink_glm_job(geno_dir,geno_pre,pheno_file,covar_file,out_dir,out_pre)


if (__name__=="__main__"):
     main()
