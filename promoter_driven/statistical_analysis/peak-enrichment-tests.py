#! /user/bin/python -O
# Jason Matthew Torres
'''
Python script for optimizing peak parameters based on credible set enrichment
Note that 256 is the average size of a DpnII restriction fragment
Usage:
module add bedtools
python script fname out_file
'''

import os,sys
import subprocess as sp
import numpy

stat_dir = "/t1-data/user/hugheslab/jtorres/analysis/wtsa/promoter_driven/statistical_analysis/"
out_dir = stat_dir + "output_files/"
ref_dir = "/t1-data/user/hugheslab/jtorres/analysis/wtsa/promoter_driven/reference_files/"
genome_file = ref_dir + "hg19.chrom.sizes"

# islet.TFs.incl_chr.bed
# ENGAGE_FG_1000G.bed
# ENGAGE_FIadjBMI1000G.bed

def get_intersect(gwas_bed,capc_bed,temp_name):
	command = ["bedtools", "intersect", "-wa","-a",gwas_bed,"-b",capc_bed,"|","uniq",">",out_dir+"inter.temp."+temp_name+".bed"]
	sp.check_call(" ".join(command),shell=True)
	command = ["cat", out_dir+"inter.temp."+temp_name+".bed","|","wc -l"]
	p = sp.Popen(" ".join(command),shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
	output, err = p.communicate(b"input data that is passed to subprocess' stdin")
	count = output.strip()
	os.remove(out_dir+"inter.temp."+temp_name+".bed")
	return count

def enrich(fname,gwas_bed):
	temp_name = "enrichment-tests"
	fin = open(fname,'r')
	fin.readline() # header
	fout = open(out_dir+"temp1."+temp_name+".bed",'w')
	for line in fin:
		l = line.strip().split()
		chrom,start,end,name =l[1],l[2],l[3],l[0]
		fout.write("\t".join([chrom,start,end,name])+"\n")
	fin.close()
	fout.close()
	command1 = ["bedtools", "sort -i", out_dir+"temp1."+temp_name+".bed", ">", out_dir+"temp2."+temp_name+".bed"]
	command2 = ["bedtools", "merge -i", out_dir+"temp2."+temp_name+".bed",">", out_dir+"temp."+temp_name+".bed"]
	sp.check_call(" ".join(command1),shell=True)
	sp.check_call(" ".join(command2),shell=True)
	os.remove(out_dir+"temp1."+temp_name+".bed")
	os.remove(out_dir+"temp2."+temp_name+".bed")
	observed = int(get_intersect(gwas_bed,out_dir+"temp."+temp_name+".bed",temp_name))
	iter_list = []
	for i in range(1,1001):
		sys.stdout.write("\r%d"%i)
		sys.stdout.flush()
		command = ["bedtools", "shuffle -i", out_dir+"temp."+temp_name+".bed", "-g",genome_file, "-chrom", ">", out_dir+"shuffle.temp."+temp_name+".bed"]
		sp.check_call(" ".join(command),shell=True)
		rand_count = get_intersect(gwas_bed,out_dir+"shuffle.temp."+temp_name+".bed",temp_name)
		iter_list.append(int(rand_count))
		os.remove(out_dir+"shuffle.temp."+temp_name+".bed")
	print ("\nObserved: %d" % observed)

	pval =  (sum([x>=observed for x in iter_list])+1) / float(1000+1)
	print ("Pvalue: %f" % pval)
	return [observed,numpy.mean(iter_list), observed/numpy.mean(iter_list), pval]

def main():
	file_list = ["peakC-modeled_fdr0.05_win11.txt","peakC-modeled_fdr0.05_win42.txt","peakC-comparative_abscutoff-170_win11.txt","peakC-comparative_abscutoff-170_win42.txt"]
	afile_list = ["islet.TFs.incl_chr.bed", "FOXA2.bed","NKX2.2.bed","NKX6.1.bed","MAFB.bed","PDX1.bed","ENGAGE_FG_1000G.bed","ENGAGE_FIadjBMI1000G.bed"]
	gwas_list = os.listdir(ref_dir+"gwas_beds/")
	gwas_list = [x for x in gwas_list if "diabetes" in x]
	afile_list =  gwas_list + afile_list
	print afile_list
	out_file = stat_dir + "enrichment-testsV3.txt"
	fout = open(out_file,'w')
	fout.write("\t".join(["test","feature","observed.count","null.mean","enrich.factor","p.val"])+"\n")
	for f in file_list:
		tname = f.split(".txt")[0]
		print tname
		for a in afile_list:
			try:
				aname = a.split(".bed")[0]
				print(aname)
				if a in gwas_list:
					out_list = enrich(out_dir+f,ref_dir+"gwas_beds/"+a)
				else:
					out_list = enrich(out_dir+f,ref_dir+a)
				out_list = [str(x) for x in out_list]
				write_list = [tname,aname] + out_list
				fout.write("\t".join(write_list)+"\n")
			except:
				print aname
	fout.close()

if (__name__=="__main__"): main()
