#!/usr/bin/python -O
# Jason Matthew Torres
'''
Python wrapper for evaluating eQTL genotype information in endoB and LCL samples
Usage: python 05_assess-genotypes.py
'''

import sys, os, gzip
import subprocess as sp


vcftools="/apps/well/vcftools/0.1.14-gcc4.7.2/bin/vcftools"

eqtl_dir = "/well/got2d/jason/reference/islet/eqtls/inspire/"
eqtl_res_dir = eqtl_dir + "nominal_pass/output/"
eqtlvcf_dir = eqtl_dir+"vcf/" # inspire.vcf.gz
# /apps/well/vcftools/0.1.14-gcc4.7.2/bin/vcftools --gzvcf inspire.vcf.gz --snp  rs3107975 --recode --out test
root_dir = "/well/got2d/jason/projects/wtsa/promoter_driven/"
stat_dir= root_dir + "statistical_analysis/"
ref_dir= root_dir + "reference_files/"

job_dir = stat_dir + "jobs/"
log_dir = stat_dir + "logs/"

def write_snp_file():
	fin = open(ref_dir + "eqtl-index.bed",'r')
	fout = open(ref_dir + "eqtl-index.snps",'w')
	fin.readline()
	for line in fin:
		l = line.strip().split()
		name = l[-1]
		rsid = name.split("_")[0]
		fout.write(rsid+"\n")
	fin.close()
	fout.close()

def run_job():
	command_list = [vcftools, "--gzvcf",eqtlvcf_dir+"inspire.vcf.gz", "--snps",
					ref_dir+"eqtl-index.snps","--recode","--out", stat_dir+"genotypes/eqtl-inspire.vcf"]
	command = " ".join(command_list)

	script = '''
#$ -N vcf_subset
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %svcf_subset.error
#$ -o %svcf_subset.out
#$ -V

module load python/2.7.11
%s

	''' % (log_dir,log_dir,command)
	print ("Writing job script")
	fout = open(job_dir+"vcf_subset_job.sh",'w')
	fout.write(script)
	fout.close()
	call = "qsub " + job_dir+"vcf_subset_job.sh"
	sp.check_call(call,shell=True)

def write_snp_info_file():
	fin = open(stat_dir+"genotypes/eqtl-inspire.vcf.recode.vcf",'r')
	fin.readline()
	fout = open(stat_dir+"genotypes/eqtl-alleles.txt",'w')
	fout.write("\t".join(["CHROM","POS","ID","REF","ALT"])+"\n")
	for line in fin:
		l = line.strip().split()
		write_list = [l[0],l[1],l[2],l[3],l[4]]
		print write_list
		fout.write("\t".join(write_list)+"\n")
	fin.close()
	fout.close()

def build_summary_file():
	fin = open(stat_dir+"genotypes/"+"eqtl-alleles.txt",'r')
	fin.readline()
	a_dic = {}
	for line in fin:
		l = line.strip().split()
		rsid = l[2]
		a_dic[rsid] = l
	fin.close()
	print a_dic
	snp_list = a_dic.keys()
	print snp_list
	e_dic = {}
	for snp in snp_list:
		e_dic[snp] = []
	print e_dic
	fin = gzip.open(eqtl_res_dir+"nominal.all.chunks.txt.gz",'rb')
	#fin = gzip.open(eqtl_res_dir+"eqtls_fdr05.txt.gz",'rb')

	count = 0
	for line in fin:
		count += 1
		sys.stdout.write("\r%d"%count)
		sys.stdout.flush()
		l = line.strip().split()
		rsid = l[1]
		try:
			e_dic[rsid].append(l)
		except:
			pass
	fin.close()
	print e_dic
	fout = open(stat_dir+"genotypes/"+"eqtl-summary-temp.txt",'w')
	head_list = ["CHROM","POS","ID","REF","ALT", "Gene", "rsid", "Distance", "P","Slope"]#,"Q"]
	fout.write("\t".join(head_list)+"\n")
	for snp in snp_list:
		ll = e_dic[snp]
		for i in range(0,len(ll)):
			write_list = a_dic[snp] + e_dic[snp][i]
			fout.write("\t".join(write_list)+"\n")
	fout.close()



def main():
	#write_snp_file()
	#run_job()
	#write_snp_info_file()
	build_summary_file()

if (__name__=="__main__"): main()
