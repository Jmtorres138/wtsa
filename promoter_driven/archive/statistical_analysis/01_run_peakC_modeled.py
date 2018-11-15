#! /user/bin/python -O
# Jason Matthew Torres
'''
Python script for submitting multiple peakC runs
Usage: python 01_run_peakC_modeled.py
'''
import os,sys
import subprocess as sp

stat_dir = "/t1-data/user/hugheslab/jtorres/analysis/wtsa/promoter_driven/statistical_analysis/"
out_dir = stat_dir +"output_files/"
ref_dir = "/t1-data/user/hugheslab/jtorres/analysis/wtsa/promoter_driven/reference_files/"

gwas_bed = ref_dir + "ENGAGE_FG_1000G_Bonferroni.bed"#"Islets_Active_enhancers.chromatinStates.bed" #"genetic_credible-sets_independent_ppa10.bed" #  #
write_file = out_dir + "enrich_peakC-modeled_engage-fg-1000G.txt" #"enrich_peakC-modeled_islet-enhancers.txt" #"enrich_peakC-modeled_gwasPPA-10_full.txt" #  #

def run_job(win, fdr):
	output_file = out_dir+"peakC-modeled_fdr"+str(fdr)+"_win"+str(win)+".txt"
	command = ["Rscript", "--vanilla", stat_dir+"01_peakC_modeled.R", str(win), str(fdr)]
	com = " ".join(command)
	command2 = ["python",stat_dir+"01_optimize_peakC-modeled.py",output_file,gwas_bed,write_file]
	com2 = " ".join(command2)

	job_name = str(win)+"_"+str(fdr)
	fout = open(stat_dir+"jobs/"+"job_optim_mod-"+job_name+".sh",'w')
	script = '''
#!/bin/bash -l
#$ -cwd
#$ -M jtorres
#$ -m eas
#$ -j n
#$ -N optim_mod_%s

module add bedtools

%s
%s

	''' % (job_name,com,com2)
	fout.write(script)
	fout.close()
	command = ["qsub", stat_dir+"jobs/"+"job_optim_mod-"+job_name+".sh"]
	if os.path.exists(output_file)==False:
		sp.check_call(" ".join(command),shell=True)
		#pass


def main():
	fdr_list = ["0.01","0.05","0.10","0.20"]
	for f in fdr_list:
		for w in range(1,100):
			run_job(w,f)

if (__name__=="__main__"): main()
