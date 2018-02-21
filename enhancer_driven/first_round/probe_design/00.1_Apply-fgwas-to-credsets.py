#!/usr/bin/python -O
# Jason Matthew Torres
'''
Previously ran fgwas assessing 24 islet and genome annotations in the
diagram imputed to 1K Genomes dataset, the best model had 14 annotations with
penalty = 0.5. Will now apply this model to the MetaboChip and Diagram credible
sets with fGWAS
Usage: python 00.1_Apply-fgwas-to-credsets.py
'''
# libraries
import sys,os,gzip
import subprocess as sp
import operator
import time
from select import select

# globals
fgwas = "/users/mccarthy/jmtorres/software/fgwas-0.3.6/bin/fgwas"
model_list = ["islet_state8","cds","islet_state5","islet_atac","islet_state2"]
ridgeparam = str(0.65)
cur_dir = "/well/got2d/jason/projects/wtsa/enhancer_driven/probe_design/"
output_dir = "/well/got2d/jason/projects/wtsa/enhancer_driven/probe_design/fgwas_output/"

# functions

def apply_to_metabochip():
    home_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/credsets_metabochip/fgwas_input/"
    input_file = home_dir + "fgwas-metabo_input27_var100.txt.gz"
    #input_file = home_dir + "fgwas-metabo_input27.txt.gz"
    keep_list = list(model_list)
    job_file = cur_dir + "jobs/job_metabo.sh"
    fout=open(job_file,'w')
    if "distance_tss" in model_list:
        keep_list.remove("distance_tss")
        model_sub = "+".join(keep_list)
        command_list = [fgwas, "-i", input_file, "-cc", "-fine",
                        "-dists", "distance_tss"+":"+home_dir+"dist_model",
                        "-w", model_sub,
                        "-p", ridgeparam, "-xv", "-print",
                        #"-o", output_dir+"metabochip_gwasmod"]
                        "-o", output_dir+"metabochip_gwasmod_var100"]

    else:
        command_list = [fgwas, "-i", input_file, "-cc", "-fine",
                        "-w", "+".join(keep_list), "-p", ridgeparam, "-xv", "-print",
                        #"-o", output_dir+"metabochip_gwasmod"]
                        "-o", output_dir+"metabochip_gwasmod_var100"]

    command = " ".join(command_list)
    script='''
#$ -N job_metab
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
#$ -V

%s
    ''' % (cur_dir+"logs/","job_metab",
    cur_dir+"logs/","job_metab", command)
    fout.write(script)
    fout.close()
    call = ["qsub", job_file]
    sp.check_call(call)

def apply_to_diagram():
    home_dir = "/well/got2d/jason/projects/t2d-integration/fgwas/credsets_diagram_1KG/fgwas_input/"
    input_file = home_dir + "fgwas_input27_var100.txt.gz"
    keep_list = list(model_list)
    job_file = cur_dir + "jobs/job_diag.sh"
    fout=open(job_file,'w')
    if "distance_tss" in model_list:
        keep_list.remove("distance_tss")
        model_sub = "+".join(keep_list)
        command_list = [fgwas, "-i", input_file, "-cc", "-fine",
                        "-dists", "distance_tss"+":"+home_dir+"dist_model",
                        "-w", model_sub,
                        "-p", ridgeparam, "-xv", "-print",
                        #"-o", output_dir+"diagram_gwasmod"]
                        "-o", output_dir+"diagram_gwasmod_var100"]

    else:
        command_list = [fgwas, "-i", input_file, "-cc", "-fine",
                        "-w", "+".join(keep_list), "-p", ridgeparam, "-xv", "-print",
                        #"-o", output_dir+"diagram_gwasmod"]
                        "-o", output_dir+"diagram_gwasmod_var100"]

    command = " ".join(command_list)
    script='''
#$ -N job_diag
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s%s.error
#$ -o %s%s.out
#$ -V

%s
    ''' % (cur_dir+"logs/","job_diag",
    cur_dir+"logs/","job_diag", command)
    fout.write(script)
    fout.close()
    call = ["qsub", job_file]
    sp.check_call(call)



def main():
    apply_to_metabochip()
    apply_to_diagram()

if (__name__=="__main__"): main()
