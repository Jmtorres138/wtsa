#! /user/bin/python -O
# Jason Matthew Torres
'''
Python script for extracting relevant summary gwas input from Mahajan et al. \
    2018b (i.e. conditioned summary statistics), that will be needed in \
    subsequent colocalisation analyses with COLOC (Chris Wallace)
Usage:
module load Python/3.7.4-GCCcore-8.3.0
python 5.1_prepare-gwas-input-coloc.py 
'''

#import os,sys,gzip
#import subprocess as sp

import pandas as pd 
import numpy as np 


serv_dir1 = "/well/got2d/jason/"#"/home/jason/science/servers/FUSE/"
serv_dir2 = "/well/mccarthy/users/jason/"#"/home/jason/science/servers/FUSE5/"

work_dir = serv_dir2+"projects/wtsa/joint_analyses/"
out_dir = work_dir + "coloc_analysis_files/conditioned_gwas_summarystats/"
cond_dir = serv_dir1+"projects/t2d-integration/fgwas/diagram_hrc/"+\
"ukbb-diamante-euro-manuscript/conditional/fgwas_input_files/"

cred_file = serv_dir2 + "datasets/diamante_hrc/"+\
    "genetic_credible_sets/gencred.txt"
    
cred_df = pd.read_csv(cred_file,sep="\t")

id_list = sorted(list(set(cred_df['CondID'])))
special_list = ['132_1','132_2','132_3','132_4','132_5',
                '133_1','133_2','133_3','133_4','133_5',
                '133_6','133_7','133_8','133_9','133_10',
                '210_1','210_2','86_1','87_1']



# Read in fgwas data frames 

cond1_df = pd.read_csv(cond_dir+"ukbb_diamante-euro.cond1.fgwas.gz",\
                       sep="\t")
cond2_df = pd.read_csv(cond_dir+"ukbb_diamante-euro.cond2.fgwas.gz",\
                       sep="\t")
cond3_df = pd.read_csv(cond_dir+"ukbb_diamante-euro.cond3.fgwas.gz",\
                       sep="\t")
cond4_df = pd.read_csv(cond_dir+"ukbb_diamante-euro.cond4.fgwas.gz",\
                       sep="\t")
cond5_df = pd.read_csv(cond_dir+"ukbb_diamante-euro.cond5.fgwas.gz",\
                       sep="\t")    
cond6_df = pd.read_csv(cond_dir+"ukbb_diamante-euro.cond6.fgwas.gz",\
                       sep="\t")
cond7_df = pd.read_csv(cond_dir+"ukbb_diamante-euro.cond7.fgwas.gz",\
                       sep="\t")
cond8_df = pd.read_csv(cond_dir+"ukbb_diamante-euro.cond8.fgwas.gz",\
                       sep="\t")

    
    
    
def extract_signal_info(sig,range=500000):
    is_sig = cred_df['CondID']==sig
    sub_df = cred_df[is_sig]
    lead_rsid = str(list(set(sub_df['lead.rsid']))[0])
    index_snp = str(list(set(sub_df['IndexSNP']))[0])
    chrom,pos = index_snp.split(":")[0],int(index_snp.split(":")[1])
    win_start = pos - range
    win_end = pos + range
    symbol = str(list(set(sub_df['symbol']))[0])
    return([sig,symbol,lead_rsid,chrom,pos,win_start,win_end])


def collect_summary_stats(sig):
    l = extract_signal_info(sig)
    cond_num = sig.split("_")[1]
    if sig in special_list:
        #print("Special")
        file_path = cond_dir + "ukbb_diamante-euro."+sig+".fgwas.gz"
        data_df = pd.read_csv(file_path,sep="\t")
    elif cond_num=='1':
        data_df=cond1_df
    elif cond_num=='2':
        data_df=cond2_df
    elif cond_num=='3':
        data_df=cond3_df
    elif cond_num=='4':
        data_df=cond4_df  
    elif cond_num=='5':
        data_df=cond5_df    
    elif cond_num=='6':
        data_df=cond6_df
    elif cond_num=='7':
        data_df=cond7_df
    elif cond_num=='8':
        data_df=cond8_df
    else:
        print("Conditional signal number is not present")
    a = np.array(data_df['CHR']==l[3]) # chromosome must match 
    b = np.array(data_df['POS']>=l[-2]) # must be within range of index snp 
    c = np.array(data_df['POS']<=l[-3]) # must be within range of index snp
    keep_list =  a & b & c 
    sub_df = data_df[keep_list][['CHR','POS', 'SNPID', 'F', 'Z', 'PVAL', 
                                 'NCASE', 'NCONTROL']]
    write_name = out_dir + "signal-"+sig+"-"+l[1]+"-"+l[2]+"-"+l[3]+"-"+\
        str(l[4])+".csv"
    sub_df.to_csv(write_name,index=False)
    

def main():
    for sig in id_list:
        print(sig)
        collect_summary_stats(sig)

if (__name__=="__main__"): 
    main()
