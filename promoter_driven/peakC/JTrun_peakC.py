#! /user/bin/python -O 
# Jason Matthew Torres
'''
Python script for submitting multiple peakC runs 
Usage: python JTrun_peakC.py 
'''

import subprocess as sp 

peakc_dir = "/home/wellgen/jtorres/projects/promoter-driven/peakC_interactions/"


def run_job(win, fdr):
	command = ["Rscript", "--vanilla", peakc_dir+"run_peakC.R", str(win), str(fdr)]
	com = " ".join(command)
	print com
	sp.check_call(com,shell=True)

def main():
	fdr_list = ["0.01"]#,"0.20"]#,"0.05","0.10","0.20"]
	for f in fdr_list:
		for w in range(11,12):
			run_job(w,f)

if (__name__=="__main__"): main() 