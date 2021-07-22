

pipePath="/home/molhaem2/telenius/CCseqBasic/CS5"
Sample="Endo_A"
Genome="hg19"
OligoFile="/t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/second_round/rerun_pipeline/OligoFileRedo.txt"
PublicPath="/public/jtorres/wtsa/enhancer_driven/second_round/rerun_pipeline/${Sample}"


${pipePath}/pipe.sh -o "${OligoFile}" -s "${Sample}" --pf "${PublicPath}" --genome ${Genome} --R1 /t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/second_round/fastq_files/original_sequencing/Endo_A_R1.fastq --R2 /t1-data/user/hugheslab/jtorres/wtsa/enhancer_driven/second_round/fastq_files/original_sequencing/Endo_A_R2.fastq --onlyBlat
