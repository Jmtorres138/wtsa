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

stat_dir = "/t1-data/user/hugheslab/jtorres/analysis/wtsa/promoter_driven/statistical_analysis/"
out_dir = stat_dir + "output_files/"
ref_dir = "/t1-data/user/hugheslab/jtorres/analysis/wtsa/promoter_driven/reference_files/"
genome_file = ref_dir + "hg19.chrom.sizes"

fname = sys.argv[1] # fname = out_dir + "peakC-modeled_fdr0.01_win11.txt"
gwas_bed = sys.argv[2] # ref_dir + "genetic_credible-sets_independent_ppa01.bed"
out_file = sys.argv[3] # out_dir + "enrich_peakC-modeled.txt"

def get_intersect(gwas_bed,capc_bed,temp_name):
	command = ["bedtools", "intersect", "-wa","-a",gwas_bed,"-b",capc_bed,"|","uniq",">",out_dir+"inter.temp."+temp_name+".bed"]
	sp.check_call(" ".join(command),shell=True)
	command = ["cat", out_dir+"inter.temp."+temp_name+".bed","|","wc -l"]
	p = sp.Popen(" ".join(command),shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
	output, err = p.communicate(b"input data that is passed to subprocess' stdin")
	count = output.strip()
	os.remove(out_dir+"inter.temp."+temp_name+".bed")
	return count

def enrich(fname,fdr,win):
	temp_name = "fdr"+str(fdr)+"_win"+str(win)
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
	#command2 = ["bedtools", "merge -i", out_dir+"temp2.bed", "-d 256",">", out_dir+"temp3.bed"]
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
	return [observed,pval]

def main():
	if os.path.exists(out_file)==False:
		fout = open(out_file,'w')
		fout.write("\t".join(["fdr","window","gwas.count","p.val"])+"\n")
	else:
		fout = open(out_file,'a')
	fdr = fname.split("fdr")[1].split("_win")[0]
	win = fname.split("win")[1].split(".txt")[0]
	gcount = enrich(fname,fdr,win)[0]
	pval = enrich(fname,fdr,win)[1]
	fout.write("\t".join([fdr,win,str(gcount),str(pval)])+"\n")
	fout.close()


if (__name__=="__main__"): main()
