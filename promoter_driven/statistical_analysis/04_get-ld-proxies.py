#!/usr/bin/python -O
# Jason Matthew Torres
'''
Python wrapper for running ld analysis
Usage: python JT_ld-analysis.py
'''

import sys, os, gzip
import subprocess as sp


vcftools="/apps/well/vcftools/0.1.14-gcc4.7.2/bin/vcftools"
ref_dir="/well/got2d/jason/reference/1KGenomes/"
root_dir = "/well/got2d/jason/projects/wtsa/promoter_driven/"
stat_dir= root_dir + "statistical_analysis/"
promref_dir= root_dir + "reference_files/"

job_dir = stat_dir + "jobs/"
log_dir = stat_dir + "logs/"
ld_dir = stat_dir + "ld-proxies/"


def make_pos_files():
	sys.stdout.write("Creating position files...\n")
	fin = open(promref_dir+"eqtl-index.bed",'r')
	head_list = ["#chr","pos"]
	for line in fin:
		l = line.strip().split()
		chrom,pos,name = l[0].split("chr")[1],l[2],l[3]
		fout = open(ld_dir+name+".positions",'w')
		fout.write("\t".join(head_list)+"\n")
		fout.write("\t".join([chrom,pos])+"\n")
		fout.close()
	fin.close()


def subset_ld_results(thresh=0.80):
	hapfiles = [x for x in os.listdir(ld_dir) if ".list.hap.ld" in x]
	for f in hapfiles:
		name = f.split(".list.hap.ld")[0]
		print ("\n" + name)
		fin = open(ld_dir+f,'r')
		fout = open(ld_dir+name+".ld-proxies.txt",'w')
		head_list = fin.readline().strip().split()
		fout.write("\t".join(head_list)+"\n")
		count = 0
		for line in fin:
			count+=1
			sys.stdout.write("\r%d"%count)
			sys.stdout.flush()
			l = line.strip().split()
			r2 = float(l[-1])
			if r2 >= thresh:
				fout.write("\t".join(l)+"\n")
		fin.close()
		fout.close()
	print("\n" + "Process Complete")


def run_job(position_list_file):
	name = position_list_file.split("/")[-1].split(".positions")[0]
	fin = open(position_list_file,'r')
	fin.readline()
	l = fin.readline().strip().split()
	chrom = "chr" + str(l[0])
	pos = l[1]
	chromosome = chrom[3:]
	mystart = int(pos) - 1e6
	myend = int(pos) + 1e6
	command_list = [vcftools, "--gzvcf", ref_dir+"ALL."+str(chrom)+".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
	"--chr", str(chromosome), "--from-bp", str(mystart), "--to-bp", str(myend), "--remove-indels",
	"--hap-r2-positions", position_list_file,  "--ld-window-bp", "1000000",
	"--out", ld_dir+name]
	command = " ".join(command_list)

	script = '''
#$ -N %s_vcf_ld_proxies
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s_vcf_ld.error
#$ -o %s%s_vcf_ld.out
#$ -V

module load python/2.7.11
%s

	''' % (name,
		   log_dir,name,
		   log_dir,name,
		   command)
	print ("Writing job script")
	fout = open(job_dir+name+"_job.sh",'w')
	fout.write(script)
	fout.close()
	call = "qsub " + job_dir+name+"_job.sh"
	sp.check_call(call,shell=True)

def main():
	#make_pos_files()
	#posfiles = [x for x in os.listdir(ld_dir) if ".positions" in x]
	#for f in posfiles:
	#	fpath = ld_dir + f
	#	print f
	#	run_job(fpath)
	subset_ld_results()

if (__name__=="__main__"): main()
