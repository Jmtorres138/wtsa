import sys,os

cur_dir = "/t1-data/user/hugheslab/jtorres/wtsa/promoter_driven/"

samp_list = ["BlymphA","BlymphB","BlymphC","EndoBA","EndoBC","EndoBD"]


for samp in samp_list:
    samp_dir = cur_dir + samp + "/"
    if os.path.exists(samp_dir) == False:
        os.makedirs(samp_dir)
        fout = open(samp_dir + "PIPE_fastqPaths.txt",'w')
        count = len(os.listdir(cur_dir+"fastq_files/" + samp+"_R1"))
        for i in range(0,count):
            a = cur_dir +"fastq_files/" + samp + "_R1/" + samp + "_R1_" + str(i) + ".fastq.gz"
            b = cur_dir +"fastq_files/" + samp + "_R2/" + samp + "_R2_" + str(i) + ".fastq.gz"
            if os.path.isfile(a) == True and os.path.isfile(b)==True:
                fout.write(" ".join([a,b])+"\n")
            else:
                print a#"Problem with sample %s : count %i" % (samp,i)
        fout.close()
    else:
        print("Directory already exists: Must delete to create new fastqPaths files")
