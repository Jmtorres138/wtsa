#!/usr/bin/python -O
# Jason Matthew Torres
'''
Extract read counts
Usage: python 00.1_script
'''
# libraries
import sys,os
import gzip

samp_name = sys.argv[1] # "Endo_A"
cur_dir = "/well/mccarthy/users/jason/projects/wtsa/enhancer_driven/first_round/evaluate_pipeline_output/"
report_dir = "/well/mccarthy/projects/ng_capture-C/endoC-BH1/enhancer_driven/first_round/bioinformatics_pipeline/" + samp_name + "/F6_greenGraphs_combined_"+samp_name+"_CS5/"
report_file = report_dir + "COMBINED_report_CS5.txt"

def build_result_dic(rep_file,samplename):
    dic = {}
    fin = open(rep_file,'r')
    for line in fin:
        l = line.strip().split()
        try:
            if l[1] == "12": # Reporters before final filtering steps
                #print line
                dic[l[0]] = [l[-1]]
            if l[1] == "15": # Capture fragments (final count)
                dic[l[0]].append(l[-1])
            if l[1] == "17a": # Reporter fragments (final count)
                dic[l[0]].append(l[-1])
            if l[1] == "17b": # Reporter fragments CIS (final count)
                dic[l[0]].append(l[-1])
                dic[l[0]].append(samplename)
        except:
            pass
    fin.close()
    return dic

def dic2file(dic,outname):
    fout = open(outname,'w')
    head_list = ["capture.name","all.reporters","all.captures","all.reporters.filtered","cis.reporters.filtered","sample"]
    fout.write("\t".join(head_list)+"\n")
    for key in dic:
        write_list = [key] + dic[key]
        fout.write("\t".join(write_list)+"\n")
    fout.close()

def main():
    endo_dic = build_result_dic(report_file,samp_name)
    dic2file(endo_dic,cur_dir+"report-"+samp_name+".txt")

if (__name__=="__main__"): main()
