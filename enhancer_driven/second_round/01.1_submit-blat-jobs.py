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
work_dir = serv_dir+"projects/wtsa/enhancer_driven/second_round/"
blat_dir = work_dir + "blat/"
log_dir = blat_dir + "logs/" #
job_dir = blat_dir + "jobs/" #
oligo_dir = blat_dir + "oligo_files/"
oligo_file = work_dir + "input_files/Oligo-file_merged_GRCh38_no-chr.txt"
fastq_dir = work_dir + "fastq_files/original_sequencing/"
fastq_prefix = "Endo_A_R" # example of full name expected: Endo_A_R1.fastq.gz
out_dir = blat_dir + "output_files/"
reuse_dir = work_dir + "01_REUSE_blat/"

def write_blat_script(cap):
    script = '''#!/bin/bash
pipePath="/well/mccarthy/users/jason/projects/wtsa/software/CCseqBasicS-VS1.0.8/"
pwd
Read1="%s1.fastq.gz"
Read2="%s2.fastq.gz"
Genome="hg38"
OligoFile="%s"
PublicPath="%spublic/"
Sample="EndoB_rep1"
Blatpath="%s01_REUSE_blat/"
rundir=$( pwd )
echo "Running pipe.sh in folder ${rundir}"
echo "${pipePath}/CCseqBasic5.sh.sh -o "${OligoFile}" -s "${Sample}" --pf "${PublicPath}" --genome ${Genome} --chunkmb 1012 --R1 ${Read1} --R2 ${Read2}  --gz --BLATforREUSEfolderPath ${Blatpath} --onlyBlat"
${pipePath}/CCseqBasic5.sh -o "${OligoFile}" -s "${Sample}" --pf "${PublicPath}" --genome ${Genome} --chunkmb 1012 --R1 ${Read1} --R2 ${Read2} --gz  --BLATforREUSEfolderPath ${Blatpath} --onlyBlat
echo "All done !"
echo
    ''' % (fastq_dir+fastq_prefix,fastq_dir+fastq_prefix,oligo_dir+cap+".txt",
    work_dir,work_dir)
    fout = open(job_dir+cap+".blat.sh",'w')
    fout.write(script)
    fout.close()

def submit_job(cap):
    command0 = ["mkdir",out_dir+cap]
    if os.path.isdir(out_dir+cap)==False:
        sp.check_call(command0)
    script = '''#!/bin/bash
#$ -P mccarthy.prjc
#$ -N blat-%s
#$ -q short.qc
#$ -pe shmem 1
#$ -o %s.out
#$ -e %s.err
#$ -wd %s

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

sh %s
    ''' % (cap,log_dir+cap,log_dir+cap,out_dir+cap,
    job_dir+cap+".blat.sh")
    fout = open(job_dir+cap+".job.sh",'w')
    fout.write(script)
    fout.close()
    command = ["qsub",job_dir+cap+".job.sh"]
    sp.check_call(command)

def submit_job_long(cap):
    command0 = ["mkdir",out_dir+cap]
    if os.path.isdir(out_dir+cap)==False:
        sp.check_call(command0)
    script = '''#!/bin/bash
#$ -P mccarthy.prjc
#$ -N blat-%s
#$ -q long.qc
#$ -pe shmem 2
#$ -o %s.out
#$ -e %s.err
#$ -wd %s

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

sh %s
    ''' % (cap,log_dir+cap,log_dir+cap,out_dir+cap,
    job_dir+cap+".blat.sh")
    fout = open(job_dir+cap+".job.sh",'w')
    fout.write(script)
    fout.close()
    command = ["qsub",job_dir+cap+".job.sh"]
    sp.check_call(command)

def run_blat_jobs():
    fin = open(oligo_file,'r')
    for line in fin:
        l = line.strip().split()
        cap = l[0]
        fout = open(oligo_dir+cap+".txt",'w')
        fout.write("\t".join(l)+"\n")
        fout.close()
        write_blat_script(cap)
        submit_job(cap)
    fin.close()

def main():
    #run_blat_jobs() # Run once then run remaining code when all jobs completed
    fin = open(oligo_file,'r')
    for line in fin:
        l = line.strip().split()
        cap = l[0]
        blat_file = reuse_dir + "TEMP_"+cap+"_blat.psl"
        if os.path.isfile(blat_file)==False:
            print(cap)
            submit_job_long(cap)
    fin.close()


if (__name__=="__main__"):
    main()
