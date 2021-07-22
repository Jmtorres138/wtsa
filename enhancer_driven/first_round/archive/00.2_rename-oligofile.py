# quick checks for character length
import re
infile = "/t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/first_round/renamed_oligos.txt"
infile2 = "/t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/first_round/Oligofile-temp.txt"
outfile = "/t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/first_round/renamed_oligos_full.txt"
outfile2 = "/t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/first_round/Oligofile.txt"

fin = open(infile,'r')
fin2 = open(infile2,'r')
fout = open(outfile,'w')
fout2 = open(outfile2,'w')
dic = {}
for line in fin:
    l = line.strip().split()
    if l[0][0] == "r":
        dic[l[0]] = l[0]
        fout.write("\t".join([l[0],l[0]])+"\n")
    else:
        s = "_".join("__".join(l[0].split("__")[0:4]).split("_")[0:8])
        dic[l[0]] = s
        fout.write("\t".join([l[0],s])+"\n")
        #print dic
fin.close()
fout.close()


for line in fin2:
    l = line.strip().split()
    name = l[0]
    try:
        new_name = dic[name]
        l[0] = new_name
        fout2.write("\t".join(l)+"\n")
    except:
        fout2.write("\t".join(l)+"\n")
fin2.close()
fout2.close()
