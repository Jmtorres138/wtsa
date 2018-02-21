#!/usr/bin/python -O
# Jason Matthew Torres
'''
Python wrapper for running ld analysis
Usage: python JTget_ld.py
'''

import sys, os, gzip
import subprocess as sp


vcftools="/apps/well/vcftools/0.1.14-gcc4.7.2/bin/vcftools"
ref_dir="/well/got2d/jason/reference/1KGenomes/"
work_dir="/well/got2d/jason/projects/wtsa/enhancer_driven/probe_design/ld/"
out_dir = work_dir + "ld_files/"
#log_dir = work_dir + "logs/"

def run_job(pos_dir, position_list_file):
	l = position_list_file.split("_")
	chrom= l[0]
	chromosome = chrom.split("chr")[1]
	pre = position_list_file.split(".positions")[0]

	command_list = [vcftools, "--gzvcf", ref_dir+"ALL."+str(chrom)+".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
	"--chr", str(chromosome), "--remove-indels", #"--min-r2", "0.01",
	"--hap-r2", "--positions", pos_dir + position_list_file,
	"--out", out_dir+pre]
	print command_list
	command = " ".join(command_list)

	script = '''
#$ -N %s_ld
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %slogs/%s_ld.error
#$ -o %slogs/%s_ld.out
#$ -V

module load python/2.7.11
%s

	''' % (pre,work_dir,pre,work_dir,pre,command)
	print ("Writing job script")
	fout = open(work_dir+"jobs/"+pre+"_Job.sh",'w')
	fout.write(script)
	fout.close()
	call = "qsub " + work_dir+"jobs/"+pre+"_Job.sh"
	sp.check_call(call,shell=True)

def main():
	posfiles = [x for x in os.listdir(work_dir+"positions/") if ".positions" in x]
	for posfile in posfiles:
			run_job(work_dir+"positions/",posfile)


if (__name__=="__main__"): main()
