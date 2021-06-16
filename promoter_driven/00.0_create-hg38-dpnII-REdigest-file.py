#! /user/bin/python -O
# Jason Matthew Torres
'''
Python script for preparing dpnII GRCh38 digest file
'''

import sys,os,re
import subprocess as sp

fasta_dir = "/well/mccarthy/users/jason/datasets/GTEx/TOPMED_reference_files/"
fasta_file = fasta_dir + "Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta"
out_dir = "/well/mccarthy/users/jason/projects/wtsa/digests/"
perl_script = "/well/mccarthy/users/jason/projects/wtsa/software/CCseqBasicS-VS1.0.8/bin/runscripts/dpnIIcutGenome4.pl"

def reformat_fasta():
    fin = open(fasta_file,'r')
    fout = open(out_dir+"hg38.fasta",'w')
    for line in fin:
        m = re.search(">",line.strip())
        if m is None:
            fout.write(line.strip()+"\n")
        else:
            l = line.strip().split()
            fout.write(l[0]+"\n")
    fin.close()
    fout.close()

def main():
    #reformat_fasta()
    command = ["perl",perl_script,out_dir+"hg38.fasta"]
    sp.check_call(command)
    # NOTE: Need to move genome coordinate file to out_dir 

if (__name__=="__main__"):
    main()
