#! /user/bin/python -O 
# Jason Matthew Torres
'''
Python script for submitting multiple peakC runs 
Usage: python 01_run_peakC_modeled.py 
'''

import subprocess as sp 

stat_dir = "/t1-data/user/hugheslab/jtorres/analysis/wtsa/promoter_driven/statistical_analysis/"


def run_job(win, fdr):
	command = ["Rscript", "--vanilla", stat_dir+"01_run_peakC_modeled.R", str(win), str(fdr)]
	com = " ".join(command)
	print com
	sp.check_call(com,shell=True)

def main():
	fdr_list = ["0.01"]#,"0.20"]#,"0.05","0.10","0.20"]
	for f in fdr_list:
		for w in range(11,12):
			run_job(w,f)

if (__name__=="__main__"): main() 
