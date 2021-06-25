#! /user/bin/python -O
# Jason Matthew Torres
'''
Usage:
module load Python/3.7.4-GCCcore-8.3.0
python script.py
'''

import os,sys,gzip
import subprocess as sp
serv_dir = "/well/mccarthy/users/jason/"
work_dir = serv_dir+"projects/wtsa/promoter_driven/"
output_dir = work_dir + "peaky_interactions/"
log_dir = output_dir + "logs/" #
job_dir = output_dir + "jobs/" #

def read_bait_list(bts_file):
    fin = open(bts_file,'r')
    header = fin.readline()
    out_list = []
    for line in fin:
        l = line.strip().split()
        bait_id = l[0]
        out_list.append(bait_id)
    fin.close()
    out_list = list(set(out_list))
    return(out_list)

def submit_job(exeriment,exp_code,bait_id):
    script = '''#!/bin/bash
#$ -P mccarthy.prjc
#$ -N peaky-%s-%s
#$ -q short.qc
#$ -pe shmem 1
#$ -o %s.out
#$ -e %s.err

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

module load R/3.6.2-foss-2019b
Rscript %s %s %s
    ''' % (exp_code,bait_id,
    log_dir+exeriment+"."+bait_id,
    log_dir+exeriment+"."+bait_id,
    work_dir+"04.0_peaky-map-interactions.R",exeriment,bait_id)
    fout = open(job_dir+exeriment+"."+bait_id + ".job.sh",'w')
    fout.write(script)
    fout.close()
    command = ["qsub",job_dir+exeriment+"."+bait_id + ".job.sh"]
    sp.check_call(command)

def run_jobs():
    exper_name = "promoter-capture"
    exper_code = "p"
    bait_list = read_bait_list(output_dir + exper_name + ".BTS.txt")
    for bait in bait_list:
        #outfile = output_dir + exper_name +"/"+ exper_name +"."+bait+".peaky-output.txt" # first chain
        outfile = output_dir + exper_name +"/chain_2/"+ exper_name +"."+bait+".peaky-output.txt" # second chain
        if os.path.isfile(outfile)==False:
            print(bait)
            submit_job(exper_name,exper_code,bait)


def main():
    run_jobs()

if (__name__=="__main__"):
    main()
