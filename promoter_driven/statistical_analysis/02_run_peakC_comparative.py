#! /user/bin/python -O 
# Jason Matthew Torres
'''
Python script for submitting multiple peakC runs 
Usage: python 02_run_peakC_comparative.py 
'''

import subprocess as sp 

stat_dir = "/t1-data/user/hugheslab/jtorres/analysis/wtsa/promoter_driven/statistical_analysis/"


def run_job(win, cutoff):
	command = ["Rscript", "--vanilla", stat_dir+"02_peakC_comparative.R", str(win), str(cutoff)]
	com = " ".join(command)
	print com
	sp.check_call(com,shell=True)

def main():
	cutoff_list = ["100"]#,"0.20"]#,"0.05","0.10","0.20"]
	for c in cufoff_list:
		for w in range(11,12):
			run_job(w,c)

if (__name__=="__main__"): main() 
