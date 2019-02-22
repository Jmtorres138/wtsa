# quick checks for character length
import re
import os.path
from shutil import copyfile
infile = "/t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/first_round/renamed_oligos_full.txt"
blat_dir = "/t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/first_round/REUSE_blat/"
blat_dir2 = "/t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/first_round/REUSE_blat_renamed/"

fin = open(infile,'r')
dic = {}
for line in fin:
    l = line.strip().split()
    dic[l[0]] = l[1]
fin.close()

f_list = os.listdir(blat_dir)#[0:1]
for f in f_list:
    copyfile(blat_dir+f,blat_dir2+f)
    new_name = re.sub("TEMP_","",f)
    new_name = re.sub("upstream_control","uc",new_name)
    new_name = re.sub("downstream_control","dc",new_name)
    new_name = re.sub("credible_variant","cv",new_name)
    new_name = re.sub("DIAMANTE_Morris","te",new_name)
    new_name = re.sub("_blat.psl","",new_name)
    try:
        rename = dic[new_name]
        rename = "TEMP_" + rename + "_blat.psl"
        print rename
    except:
        rename = "TEMP_" + new_name + "_blat.psl"
    os.rename(blat_dir2+f,blat_dir2+rename)
