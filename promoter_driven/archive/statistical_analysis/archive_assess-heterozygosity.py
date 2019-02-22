#!/usr/bin/python -O
# Jason Matthew Torres
'''
Python wrapper for evaluating eQTL genotype information in endoB and LCL samples, run on CBRG
Usage: python 05.2_assess-heterozygosity.py
'''

import sys, os, gzip
import subprocess as sp

proc_dir = "/t1-data/user/hugheslab/jtorres/promoter-driven/Processed/"
root_dir = "/t1-data/user/hugheslab/jtorres/analysis/wtsa/promoter_driven/"
stat_dir= root_dir + "statistical_analysis/"
ref_dir= root_dir + "reference_files/"
bamtemp_dir = stat_dir + "bam_files/"
job_dir = stat_dir+"jobs/"

def get_sample_list():
	samp_list = os.listdir("/t1-data/user/hugheslab/jtorres/promoter-driven/Processed/")
	samp_list.sort()
	return(samp_list)

def get_cap_dic():
	dic = {}
	fin = open(ref_dir+"eqtl-index.bed",'r')
	for line in fin:
		l = line.strip().split()
		chrom,pos,name = l[0],l[2],l[3]
		gene = name.split("_")[1]
		dic[name] = [gene,chrom+":"+pos+"-"+pos]
	fin.close()
	return(dic)


def run_job_view(capture_gene, sample, region):
	path_to_bamfile = proc_dir+sample+"/F6_greenGraphs_combined_"+sample+"_CB4/COMBINED_CB4_"+capture_gene+".bam"
	out_file_prefix = bamtemp_dir + sample + "_" + capture_gene
	#bed_file = ref_dir + "eqtl-index.bed"
	command_list1 = ["samtools", "sort", path_to_bamfile, "-o", out_file_prefix + ".bam"]
	command_list2 = ["samtools", "index", out_file_prefix + ".bam"]
	command_list3 = ["samtools", "view", out_file_prefix + ".bam", region,
					 "-b", "-o", out_file_prefix + "_chrom.bam"]#"_eqtl.bam"]
	#command_list3 = ["samtools", "view", "-L", out_file_prefix + ".bam", "-L", bed_file,
	#				 "-b", "-o", out_file_prefix + "_chrom.bam"]
	command_list4 = ["samtools", "sort", out_file_prefix + "_chrom.bam",
					 "-o", out_file_prefix + "_chrom_sorted.bam"]
	command_list5 = ["samtools", "index", out_file_prefix + "_chrom_sorted.bam"]
	command1 = " ".join(command_list1)
	command2 = " ".join(command_list2)
	command3 = " ".join(command_list3)
	command4 = " ".join(command_list4)
	command5 = " ".join(command_list5)

	script = '''
#!/bin/bash -l
#$ -cwd
#$ -M jtorres
#$ -m eas
#$ -j n
#$ -N job_%s_%s
#$ -e %s.error
#$ -o %s.out

%s
%s
%s
%s
%s

	''' % (sample, capture_gene, job_dir+sample+"_"+capture_gene, job_dir+sample+"_"+capture_gene,command1, command2, command3,command4,command5)
	print ("Writing job script")
	fout = open(job_dir+sample+"_"+capture_gene+"_job.sh",'w')
	fout.write(script)
	fout.close()
	call = "qsub " + job_dir+sample+"_"+capture_gene+"_job.sh"
	sp.check_call(call,shell=True)


def run_job_mpileup(capture_gene,samp_list,name,pos):
	out_file_prefix = bamtemp_dir + name

	file_list = []
	for samp in samp_list:
		file_list.append(bamtemp_dir+samp+"_"+capture_gene+"_chrom_sorted.bam")


	#EndoB_D_TCF7L2_eqtl_sorted.bam
	command_list1 = ["samtools", "mpileup"] + file_list + ["-o", out_file_prefix +"_temp"]
	command_list2 = ["grep", pos, out_file_prefix +"_temp"  ,">", out_file_prefix +".mpileup"]
	command_list3 = ["echo","hello"]#"rm *_temp* *.bam*"]
	command1 = " ".join(command_list1)
	command2 = " ".join(command_list2)
	command3 = " ".join(command_list3)

	script = '''
#!/bin/bash -l
#$ -cwd
#$ -M jtorres
#$ -m eas
#$ -j n
#$ -N job_mp_%s
#$ -e %s.error
#$ -o %s.out

%s
%s
%s

	''' % (capture_gene,job_dir+capture_gene,job_dir+capture_gene,command1,command2,command3)
	print ("Writing job script")
	fout = open(job_dir+"mpileup_"+capture_gene+"_job.sh",'w')
	fout.write(script)
	fout.close()
	call = "qsub " + job_dir+"mpileup_"+capture_gene+"_job.sh"
	sp.check_call(call,shell=True)


def siv():
	samp_list =  get_sample_list()
	cap_dic = get_cap_dic()

	sys.stdout.write("Sorting, indexing, viewing...\n")
	for eqtl in cap_dic.keys():
		e_list = cap_dic[eqtl]
		gene,reg =  e_list[0], e_list[1].split(":")[0]
		for samp in samp_list:
			print [gene,samp,reg]
			run_job_view(gene,samp,reg)

def dif():
	file_list = os.listdir(bamtemp_dir)
	sys.stdout.write("Deleting intermediate files...\n")
	remove_list = [x for x in file_list if "chrom" not in x]
	for f in remove_list:
		try:
			os.remove(bamtemp_dir+f)
		except:
			pass

def rpf():
	samp_list =  get_sample_list()
	cap_dic = get_cap_dic()
	sys.stdout.write("Running mpileup function...\n")
	for eqtl in cap_dic.keys():
		e_list = cap_dic[eqtl]
		gene,reg =  e_list[0], e_list[1]
		pos = reg.split(":")[1].split("-")[1]
		print pos
		run_job_mpileup(gene,samp_list,eqtl,pos)


def main():
	samp_list =  get_sample_list()
	cap_dic = get_cap_dic()
	#siv()
	#dif()
	#rpf()
	#os.system("cat " + bamtemp_dir + "*.mpileup > " + bamtemp_dir+"eqtls.mpileup")
	os.system("rm " + bamtemp_dir + "*.bam*" + " rs*" + " *temp")



if (__name__=="__main__"): main()
