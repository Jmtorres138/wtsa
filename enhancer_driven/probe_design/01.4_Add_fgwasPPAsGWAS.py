#!/usr/bin/python -O
# Jason Matthew Torres
'''
For the Credible Set files (DIAGRAM and MetaboChip), will add in the fgwas PPAs from the full GWAS DIAGRAM
imputed to 1KGenomes fgwas analysis
Usage: python 01.4_Add_fgwasPPAsGWAS.py
'''
# libraries
import sys,os,gzip


# globals
cur_dir = "/well/got2d/jason/projects/wtsa/enhancer_driven/probe_design/"
txt_dir = cur_dir + "txt/"
full_file = txt_dir + "diagram1KG-fgwas14annot.bfs.gz" # "test.txt.gz" #

# functions

def build_pos_dic(fname):
    '''
    Get chrom and pos for each credible snp
    '''
    dic = {}
    fin = gzip.open(fname,'rb')
    fin.readline() # header
    count=0
    for  line in fin:
        count+=1
        sys.stdout.write("\r%d" % count)
        sys.stdout.flush()
        l = line.strip().split()
        snpid,chrom,pos = l[0],l[1],l[2]
        dic[snpid] = [chrom,pos]
    print("\n")
    fin.close()
    fin2 = gzip.open(full_file,'rb')
    fin2.readline() # header
    count=0
    for line in fin2:
        count+=1
        sys.stdout.write("\r%d" % count)
        sys.stdout.flush()
        l = line.strip().split()
        name,ppa = l[0],l[9]
        try:
            dic[name].append(ppa)
        except:
            pass
    print("\n")
    fin2.close()
    return(dic)

def update(fname,outname):
    dic = build_pos_dic(fname)
    fin = gzip.open(fname,'rb')
    fout = gzip.open(outname,'wb')
    head_list = fin.readline().strip().split("\t")
    head_list.append("PPA.fgwas.full")
    fout.write("\t".join(head_list)+"\n")
    count=0
    for line in fin:
        count+=1
        sys.stdout.write("\r%d" % count)
        sys.stdout.flush()
        l = line.strip().split("\t")
        snpid = l[0]
        if len(dic[snpid])==3:
            l.append(dic[snpid][2])
        else:
            l.append("NA")
        fout.write("\t".join(l)+"\n")
    print("\n")
    fin.close()
    fout.close()


def main():
    update(txt_dir+"fgwas-cred-diag.txt.gz",txt_dir+"fgwas-cred-diag-updated.txt.gz")
    #update(txt_dir+"fgwas-cred-metab.txt.gz",txt_dir+"fgwas-cred-metab-updated.txt.gz")

if (__name__=="__main__"): main()
