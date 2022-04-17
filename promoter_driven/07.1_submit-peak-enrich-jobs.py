#!/usr/bin/python -O
# Jason Matthew Torres
'''
Submit jobs that each generate a feature contingency table for peaky peaks

Usage:
module load Python/3.7.4-GCCcore-8.3.0
python 07.1_submit-peak-enrich-jobs.py
'''
# libraries
import sys,os
import subprocess as sp
work_dir="/well/mccarthy/users/jason/projects/wtsa/promoter_driven/"
job_dir=work_dir+"jobs/enrichment_tests/"
log_dir=work_dir+"logs/enrichment_tests/"
data_dir="/well/mccarthy/users/jason/datasets/"
miguel_dir=data_dir+"miguel-escalada_processed_data/grch38/"
reg_dir=miguel_dir+"split_regulome/"
acc_chrom_dir=data_dir+"amp_cmdga_atlas/accessible_chromatin/"
scell_dir=data_dir+"amp_cmdga_atlas/single_cell/snATAC-seq/accessible_chromatin/"
scell_spec_dir=scell_dir+"cell_specific_peaks/"

bed_file_list=[miguel_dir+"atac_consistent_peaks.sorted.bed",\
miguel_dir+"ctcf_consistent_peaks_q001.sorted.bed",\
miguel_dir+"h3k27ac_consistent_peaks_q005.sorted.bed",\
miguel_dir+"h3k4me3_consistent_peaks_q005.sorted.bed",\
miguel_dir+"med1_consistent_peaks_q001.sorted.bed",\
miguel_dir+"smca1_consistent_peaks_q001.sorted.bed",\
reg_dir+"Active_enhancers_I.bed",reg_dir+"Active_enhancers_II.bed",\
reg_dir+"Active_enhancers_III.bed",reg_dir+"Active_promoters.bed",\
reg_dir+"Inactive_enhancers.bed",reg_dir+"Inactive_open_chromatin_regions.bed",\
reg_dir+"Strong_CTCF.bed",acc_chrom_dir+"DFF716MEQ.bed",\
scell_dir+"DFF605PWL.grch38.sorted.bed",scell_dir+"DFF048FES.grch38.sorted.bed",\
scell_dir+"DFF440XUP.grch38.sorted.bed",scell_dir+"DFF033YPK.grch38.sorted.bed",\
scell_dir+"DFF841CJP.grch38.sorted.bed",scell_dir+"DFF311FRL.grch38.sorted.bed",\
scell_dir+"DFF759JDT.grch38.sorted.bed",scell_dir+"DFF478BFI.grch38.sorted.bed",\
scell_dir+"DFF871GXZ.grch38.sorted.bed",scell_dir+"DFF933KSC.grch38.sorted.bed"]

feature_list=["islet-atac-peaks-miguel","ctcf-peaks-miguel","h3k27ac-peaks-miguel",\
"h3k4me3-peaks-miguel","med1-peaks-miguel","smca1-peaks-miguel",\
"Active_enhancers_I","Active_enhancers_II","Active_enhancers_III",\
"Active_promoters","Inactive_enhancers","Inactive_open_chromatin_regions",\
"Strong_CTCF","islet-atac-peaks-n38","acinar","stellate","ductal","gamma",\
"GCGhigh_alpha","GCGlow_alpha","INShigh_beta","INSlow_beta","SSThigh_delta","SSTlow_delta"]

bed_file_list=[scell_spec_dir+"intersection.bed",scell_spec_dir+"delta-specific.bed",
               scell_spec_dir+"beta-specific.bed",scell_spec_dir+"alpha-specific.bed",
               scell_spec_dir+"SSTlow_delta-specific.bed",scell_spec_dir+"SSThigh_delta-specific.bed",
               scell_spec_dir+"INSlow_beta-specific.bed",scell_spec_dir+"INShigh_beta-specific.bed",
               scell_spec_dir+"GCGlow_alpha-specific.bed",scell_spec_dir+"GCGhigh_alpha-specific.bed",
               scell_spec_dir+"gamma-specific.bed",scell_spec_dir+"ductal-specific.bed",
               scell_spec_dir+"stellate-specific.bed",scell_spec_dir+"acinar-specific.bed"]
feature_list=["intersection","delta-specific","beta-specific","alpha-specific",
    "SSTlow_delta-specific","SSThigh_delta-specific","INSlow_beta-specific","INShigh_beta-specific",
    "GCGlow_alpha-specific","GCGhigh_alpha-specific","gamma-specific",
    "ductal-specific","stellate-specific","acinar-specific"]



def submit_job(bed_file,feature_name):
    script = '''#!/bin/bash
#$ -cwd
#$ -P mccarthy.prjc
#$ -N ct-%s
#$ -q short.qc
#$ -pe shmem 1
#$ -o %s.out
#$ -e %s.err

echo "------------------------------------------------"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "------------------------------------------------"

module load R-bundle-Bioconductor/3.13-foss-2021a-R-4.1.0
WORKDIR=/well/mccarthy/users/jason/projects/wtsa/promoter_driven
Rscript $WORKDIR/07.0_peak-feature-enrichments.R %s %s
  ''' % (feature_name,log_dir+feature_name,log_dir+feature_name,bed_file,feature_name)
    fout = open(job_dir+feature_name+"_job.sh",'w')
    fout.write(script)
    fout.close()
    command = ["qsub",job_dir+feature_name+"_job.sh"]
    sp.check_call(command)

def main():
    for i in range(0,len(feature_list)):
        bed_file=bed_file_list[i]
        feature=feature_list[i]
        bed_name=bed_file.split("/")[-1]
        print(bed_name+" : "+feature)
        sys.stdout.write("Feature: %s\n" % feature_list[i])
        submit_job(bed_file,feature)

if (__name__=="__main__"):
     main()
