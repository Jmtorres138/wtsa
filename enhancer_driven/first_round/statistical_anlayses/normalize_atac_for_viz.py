#! /user/bin/python -O
# Jason Matthew Torres
'''
Python script for preparing normalised atac-seq peak annotation file for visualization
Usage:
module add bedtools
python script fname out_file
'''

import os,sys
import subprocess as sp
import numpy

#module load python/3.4.3 cd folder_with_bams
#mkdir bigWigs

#for I in *.bam; do
#base=basename $I .bam
#bamCoverage -b $I --binSize 100 --normalizeUsingRPKM -of bigwig -v -p 8 -o bigWigs/$base.bigWig

#done

work_dir = "/well/mccarthy/users/jason/projects/wtsa/enhancer_driven/first_round/statistical_anlayses/"
bam_dir = "/well/mccarthy/production/atac-seq/data/human_islets/full_merged_data/bams/"
elife_file = "/well/mccarthy/users/jason/projects/islet_atac/elife_samples_names_for_jason.txt"
out_dir = work_dir + "atac/"

def get_elife_list():
    fin = open(elife_file,'r')
    lst = []
    for l in fin:
        lst.append(l.strip())
    fin.close()
    return lst

def main():
    elife_list = get_elife_list()
    elife_list.append("HP1535")
    elife_list.append("HP1507_CMRL")
    print(elife_list)
    # Normalize the bams
    norm_list = []
    for f in elife_list:
        fname = bam_dir + f + ".bam"
        if os.path.isfile(fname) == True:
            norm_file = out_dir+f+".bigWig"
            norm_list.append(norm_file)
            command = "bamCoverage -b "  + fname + " --binSize 100 --normalizeUsingRPKM -of bigwig -v -p 8 -o " + norm_file # previously tried with binsize 100
            #print(command)
            sp.check_call(command,shell=True)
    # Summarize the bamss
    #multiBigwigSummary bins -b HP1507_CMRL.bigWig HP1535.bigWig -o results.npz --outRawCounts test.out
    command2 = "multiBigwigSummary bins -b " +  " ".join(norm_list) + " -o " + out_dir + "combined_results.npz --outRawCounts " + out_dir+"combined_results.txt"
    print(command2)
    sp.check_call(command2,shell=True)
if (__name__=="__main__"): main()
