#! /user/bin/python -O
# Jason Matthew Torres
'''
Usage:
module load Python/3.7.4-GCCcore-8.3.0
python script.py
'''

import os,sys,gzip
import subprocess as sp
serv_dir = "/well/mccarthy/users/jason/"
work_dir = serv_dir+"projects/wtsa/enhancer_driven/first_round/"
oligo_file = work_dir + "input_files/Oligo-file_merged_GRCh38.txt"
cur_dir = work_dir + "03_CaptureCompare/"
ccanalyser_dir = work_dir +"02_CCAnalyser/"
samp_list = ["EndoB_rep1","EndoB_rep1","EndoB_rep1", \
             "hESC_rep1","hESC_rep2","hESC_rep3"]
cc_file = cur_dir + "capture-compare-input_hg38.txt"

def capture_dict():
    dic = {}
    fin = open(oligo_file,'r')
    for line in fin:
        l = line.strip().split()
        cap = l[0]
        dic[cap]=[]
    fin.close()
    return(dic)

def append_status_info(cap_dic):
    for key in cap_dic.keys():
        #print(key)
        entry_list = []
        for samp in samp_list:
            gff_file = ccanalyser_dir+samp+"/"+"F6_greenGraphs_combined_"+ \
            samp+"_CS5/"+"COMBINED_CS5_"+key+".gff"
            status = os.path.isfile(gff_file)
            entry_list.append(status)
        cap_dic[key] = cap_dic[key] + entry_list
    return(cap_dic)

def all_same(items):
    return all(x == items[0] for x in items)

def identify_failed_captures(cap_dic):
    fout = open(cur_dir+"failed_captures.txt",'w')
    header_list = ["Capture"]+samp_list
    fout.write("\t".join(header_list)+"\n")
    cap_dic["TEST"]=[True,False,True,False,True,False]
    fail_dic = {}
    for key in cap_dic.keys():
        entry_list = cap_dic[key]
        if all_same(entry_list)==False:
            write_list = [key]+[str(e) for e in entry_list]
            fail_dic[key]=key
            fout.write("\t".join(write_list)+"\n")
    fout.close()
    return(fail_dic)

def filter_out_failures(fail_dic):
    os.rename(cc_file,cc_file+"_unfiltered")
    fin=open(cc_file+"_unfiltered",'r')
    fout=open(cc_file,'w')
    for line in fin:
        l = line.strip().split()
        cap = l[0]
        try:
            print(fail_dic[cap])
        except:
            fout.write("\t".join(l)+"\n")
    fin.close()
    fout.close()

def main():
    cap_dic = capture_dict()
    cap_dic = append_status_info(cap_dic)
    fail_dic = identify_failed_captures(cap_dic)
    filter_out_failures(fail_dic)

if (__name__=="__main__"):
    main()
