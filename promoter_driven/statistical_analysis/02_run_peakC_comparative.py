#! /user/bin/python -O
# Jason Matthew Torres
'''
Python script for submitting multiple peakC runs
Usage: python 02_run_peakC_comparative.py
'''

import sys,os
import subprocess as sp

stat_dir = "/t1-data/user/hugheslab/jtorres/analysis/wtsa/promoter_driven/statistical_analysis/"
out_dir = stat_dir +"output_files/"


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


def main():
	#cutoff_list = ["100"]#,"0.20"]#,"0.05","0.10","0.20"]
	for c in range(100,100001,100): #cutoff_list:
		for w in range(1,100):
			run_job(w,c)

if (__name__=="__main__"): main()
