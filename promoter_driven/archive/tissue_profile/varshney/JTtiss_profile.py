#!/usr/bin/python -O
# Jason Matthew Torres
'''
Profile interactions for tissue specificity with a key focus on the islet
Will use chromatin states from Varshney et al. 2016 study as a reference
Usage: module load python/2.7.11
       python JTtiss_profile
'''
# libraries
import sys,os,gzip
import subprocess as sp

# globals
bedtools = "/apps/well/bedtools/2.24.0/bedtools"

cur_dir = "/well/got2d/jason/projects/wtsa/promoter_driven/tissue_profile/varshney/"
varsh_dir = "/well/got2d/jason/reference/chromatin_segmentation/varshney_2016/chromatin_states/"
out_dir = "/well/got2d/jason/projects/wtsa/promoter_driven/tissue_profile/varshney/output/"
query_bed = sys.argv[1] #cur_dir+"stard10.query.bed" #
query_tissue = sys.argv[2]
tissue_list = [x.split(".")[0] for x in os.listdir(varsh_dir) if ".bed.gz" in x]

eval_list = ["Islets","SkeletalMuscle","Liver","Adipose","GM12878"]

def get_tiss_overlap(query_bed,tissue="Islets",maxiter=100):
    '''
    Determine states that overlap query bed in specified tissue
    '''
    pre = query_bed.split(cur_dir)[1].split(".query.bed")[0]
    bedfile = varsh_dir + tissue + ".chromatinStates.bed.gz"
    outname = out_dir + pre +"."+tissue+".intersect"
    command_list = [bedtools, "intersect", "-wb", "-a",
                    query_bed, "-b", bedfile, ">", outname]
    sys.stdout.write("Intersecting query file...\n")
    sp.check_call(" ".join(command_list),shell=True)
    sys.stdout.write("Getting overlapping states...\n")
    command_list = ["cat",outname,"|","cut","-f","8","|","uniq"]
    overstates = sp.check_output(" ".join(command_list),shell=True).split("\n")
    over_list = [x for x in overstates if x != '']
    over_dic = {}
    for state in over_list:
        over_dic[state] = 0
    #print over_dic
    sys.stdout.write("Generating null distribution...\n")
    for i in range(1,maxiter+1):
        sys.stdout.write("\r%d" % i)
        sys.stdout.flush()
        command_list = [bedtools,"shuffle","-i",bedfile,"-g",
                        cur_dir+"hg19.genome","-seed",str(i),">",outname+".temp"]
        sp.check_call(" ".join(command_list),shell=True)
        command_list = [bedtools, "intersect", "-wb", "-a",
                        query_bed, "-b", outname+".temp", ">", outname+".temp.iter"]
        sp.check_call(" ".join(command_list),shell=True)
        command_list = ["cat",outname+".temp.iter","|","cut","-f","8","|","uniq"]
        iterstates = sp.check_output(" ".join(command_list),shell=True).split("\n")
        iter_list = [x for x in iterstates if x != '']
        for state in iter_list:
            try:
                over_dic[state] = over_dic[state]+1
            except:
                pass
    for state in over_list:
        over_dic[state] = (float(over_dic[state]) + 1) / float(maxiter+1)
    output_file = out_dir+pre+"."+tissue+".enrich.states"
    fout = open(output_file,'w')
    fout.write("\t".join(["state","p.value"])+"\n")
    for items in over_dic.items():
        write_list = [items[0],str(items[1])]
        fout.write("\t".join(write_list)+"\n")
    fout.close()
    command_list = ["rm",out_dir+"*"+tissue+"*"+".temp*"]
    sp.check_call(" ".join(command_list),shell=True)


def main():
    get_tiss_overlap(query_bed,tissue=query_tissue,maxiter=1000)


if (__name__=="__main__"): main()
