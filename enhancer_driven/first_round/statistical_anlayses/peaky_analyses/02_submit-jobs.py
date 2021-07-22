#!/usr/bin/python -O
# Jason Matthew Torres
'''
Submit peaky jobs for each capture
Usage:
'''
# libraries
import sys,os,gzip
import subprocess as sp


# globals
work_dir = "/well/mccarthy/users/jason/projects/wtsa/enhancer_driven/first_round/statistical_anlayses/peaky_analyses/"
file_dir = work_dir + "analysis_files/"
job_dir = work_dir + "jobs/"
log_dir = work_dir + "logs/"
map_file = file_dir + "capture-bait-id-map.txt"
r_script = work_dir + "02_per-capture-peaks.R"



def submit_jobs():
    fin = open(map_file,'r')
    header = fin.readline()
    for line in fin:
        l = line.strip().split()
        capname =l[0]
        bait = l[1]
        filename = file_dir + capname + "_" + bait +".pky.txt"
        #print bait
        command_list = ["/apps/well/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0/bin/Rscript", "--vanilla", r_script, bait]
        command = " ".join(command_list)
        job_file = job_dir + "job_" + bait + ".sh"
        fout = open(job_file,'w')
        #module load gcc/5.4.0
        script='''
#$ -N %s
#$ -cwd -V
#$ -pe shmem 2
#$ -P mccarthy.prjc
#$ -q long.qc
#$ -e %s%s.error
#$ -o %s%s.out
echo "start time" `date`
%s
echo "end time" `date`
        ''' % ("job_"+bait, log_dir,"bait_"+bait,log_dir,"bait_"+bait,command)
        fout.write(script)
        fout.close()
        call = ["sh", job_file]
        if os.path.isfile(filename)==False:
            sp.check_call(call)
    fin.close()
    bait_file.close()

def main():
    submit_jobs()

if (__name__=="__main__"): main()
