import sys,os
cur_dir = "/t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/first_round/"
#hESC_C/F6_greenGraphs_combined_hESC_C_CM5/
samp_list = ["Endo_A","Endo_B","Endo_C","hESC_A","hESC_B","hESC_C"]


fdic = {}
for samp in samp_list:
    sdir = cur_dir + samp + "/F6_greenGraphs_combined_" + samp + "_CM5/"
    flist = os.listdir(sdir)
    for f in flist:
        try:
            fdic[f] = fdic[f] + 1
        except:
            fdic[f] = 1
remove_list = []
for key in fdic:
    if fdic[key] != 6:
        remove_list.append(key)
print len(remove_list)

for samp in samp_list:
    print samp
    sdir = cur_dir + samp + "/F6_greenGraphs_combined_" + samp + "_CM5/"
    flist = os.listdir(sdir)
    for f in flist:
        if f in remove_list:
            os.remove(sdir+f)


keep_list = []
for key in fdic:
    if fdic[key] == 6:
        keep_list.append(key)

keep_list = [s.split("COMBINED_CM5_")[1].split(".gff")[0] for s in keep_list]
fin = open(cur_dir + "capture_compare_parameters.txt",'r')
fout = open(cur_dir + "capture_compare_parameters_pruned.txt", "w")

for line in fin:
    l = line.strip().split()
    if l[0] in keep_list:
        fout.write("\t".join(l)+"\n")
fin.close()
fout.close()

print("Done")
