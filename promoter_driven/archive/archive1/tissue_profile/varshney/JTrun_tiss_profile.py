#!/usr/bin/python -O
# Jason Matthew Torres
'''
Profile interactions for tissue specificity with a key focus on the islet
Will use chromatin states from Varshney et al. 2016 study as a reference
Usage: module load python/2.7.11
       python JTtiss_profile
'''
# libraries
import sys,os,gzip
import subprocess as sp

# globals
bedtools = "/apps/well/bedtools/2.24.0/bedtools"

cur_dir = "/well/got2d/jason/projects/wtsa/promoter_driven/tissue_profile/varshney/"
varsh_dir = "/well/got2d/jason/reference/chromatin_segmentation/varshney_2016/chromatin_states/"
out_dir = "/well/got2d/jason/projects/wtsa/promoter_driven/tissue_profile/varshney/output/"
job_dir = cur_dir + "jobs/"
query_bed =  cur_dir+"stard10.query.bed" #sys.argv[1]

eval_list = ["Islets","SkeletalMuscle","Liver","Adipose","GM12878"]


def submit_job(query_bed,tissue):
    pre = query_bed.split(cur_dir)[1].split(".query.bed")[0]
    job_file = job_dir + pre+"."+tissue+".job"
    fout = open(job_file,'w')
    command_list = ["python",cur_dir+"JTtiss_profile.py",query_bed,tissue]
    command = " ".join(command_list)
    script='''
#$ -N job_diag
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
#$ -V

module load python/2.7.11
%s
    ''' % (cur_dir+"logs/","job_"+pre+"."+tissue,
    cur_dir+"logs/","job_"+pre+"."+tissue, command)
    fout.write(script)
    fout.close()
    call = ["qsub", job_file]
    #sp.check_call(call)

def main():
    for tiss in eval_list:
        submit_job(query_bed,tiss)

if (__name__=="__main__"): main()
