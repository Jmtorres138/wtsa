# quick checks for character length
import re
infile = "/t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/first_round/Oligofile-excl-cap-only.txt"#"Oligofile-original.txt"
outfile = "/t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/first_round/Oligofile-temp.txt"
outfile2 = "/t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/first_round/renamed_oligos.txt"
fin = open(infile,'r')
fout = open(outfile,'w')
fout2 = open(outfile2,'w')
for line in fin:
    l = line.strip().split()
    name = line.strip().split()[0]
    new_name = re.sub("downstream_control","dc",name)
    new_name = re.sub("upstream_control","uc",new_name)
    new_name = re.sub("credible_variant","cv",new_name)
    new_name = re.sub("DIAMANTE_Morris","te",new_name)

    l[0] = new_name
    fout.write("\t".join(l)+"\n")
    if len(new_name) >= 50:
        fout2.write(new_name+"\n")
fin.close()
fout.close()
fout2.close()
