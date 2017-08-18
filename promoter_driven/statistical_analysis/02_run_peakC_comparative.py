#! /user/bin/python -O
# Jason Matthew Torres
'''
Python script for submitting multiple peakC runs
Usage: python 02_run_peakC_comparative.py
'''

import os,sys
import subprocess as sp

stat_dir = "/t1-data/user/hugheslab/jtorres/analysis/wtsa/promoter_driven/statistical_analysis/"
out_dir = stat_dir +"output_files/"
ref_dir = "/t1-data/user/hugheslab/jtorres/analysis/wtsa/promoter_driven/reference_files/"

#gwas_bed = ref_dir + "genetic_credible-sets_independent_ppa10.bed" # "Islets_Active_enhancers.chromatinStates.bed" #
#write_file = out_dir +  "enrich_peakC-comparative_gwasPPA-10.txt" # "enrich_peakC-comparative_islet-enhancers.txt"

gwas_bed = ref_dir + "Islets_Active_enhancers.chromatinStates.bed" #
write_file = out_dir +  "enrich_peakC-comparative_islet-enhancers.txt"


def run_job(win, cutoff):


	command = ["Rscript", "--vanilla", stat_dir+"02_peakC_comparative.R", str(win), str(cutoff)]
	com = " ".join(command)
	job_name = str(win)+"_"+str(cutoff)
	fout = open(stat_dir+"jobs/"+"job_optim_comp-"+job_name+".sh",'w')
	script = '''
#!/bin/bash -l
#$ -cwd
#$ -M jtorres
#$ -m eas
#$ -j n
#$ -N optim_comp_%s

%s

	''' % (job_name,com)
	fout.write(script)
	fout.close()
	command = ["qsub", stat_dir+"jobs/"+"job_optim_comp-"+job_name+".sh"]
	if os.path.exists(out_dir+"peakC-comparative_abscutoff-"+str(cutoff)+"_win"+str(win)+".txt")==False:
		sp.check_call(" ".join(command),shell=True)


	output_file = out_dir+"peakC-comparative_abscutoff-"+str(cutoff)+"_win"+str(win)+".txt"
	command = ["Rscript", "--vanilla", stat_dir+"02_peakC_comparative.R", str(win), str(cutoff)]
	com = " ".join(command)
	command2 = ["python",stat_dir+"02_optimize_peakC-comparative.py",output_file,gwas_bed,write_file]
	com2 = " ".join(command2)

	job_name = str(win)+"_"+str(cutoff)
	fout = open(stat_dir+"jobs/"+"job_optim_comp-"+job_name+".sh",'w')
	script = '''
#!/bin/bash -l
#$ -cwd
#$ -M jtorres
#$ -m eas
#$ -j n
#$ -N optim_comp_%s

module add bedtools

%s
%s

	''' % (job_name,com,com2)
	fout.write(script)
	fout.close()
	command = ["qsub", stat_dir+"jobs/"+"job_optim_comp-"+job_name+".sh"]
	if os.path.exists(output_file)==False:
		sp.check_call(" ".join(command),shell=True)



def main():
	for c in range(10,1001,10): #cutoff_list:
		for w in [11,42,40]:
			run_job(w,c)

if (__name__=="__main__"): main()
