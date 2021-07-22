#! /user/bin/python -O
# Jason Matthew Torres
'''
module load Python/3.7.4-GCCcore-8.3.0
python 00.0_format-dpn11-reference.py
'''

import os,sys
from time import sleep
serv_dir = "/well/mccarthy/users/jason/"
dpn_file = serv_dir + "projects/wtsa/genome_dpnII_coordinates.txt"
out_file = serv_dir + "projects/wtsa/genome_dpnII_coordinates_reformatted.txt"


def reformat_file():
    count = 0
    fin = open(dpn_file,'r')
    fout = open(out_file,'w')
    for line in fin:
        count += 1
        sys.stdout.write("\r")
        sys.stdout.write("Count: %d" % count)
        sys.stdout.flush()
        l = line.strip().split(":")
        chrom = l[0]
        interval = l[1]
        l2 = interval.split("-")
        start,end = l2[0],l2[1]
        write_list = [chrom,start,end]
        fout.write("\t".join(write_list)+"\n")
    sys.stdout.write("\n")
    fin.close()
    fout.close()

def main():
    reformat_file()

if (__name__=="__main__"):
    main()
