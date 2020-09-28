"%&%" <- function(a,b) paste0(a,b)
library("data.table")
file.dir <- "/well/mccarthy/users/jason/projects/wtsa/joint_analyses/coloc_analysis_files/rasqual_output/"

args <- commandArgs(trailingOnly=TRUE)
chromo <- args[1]

print("Chromosome: " %&% chromo)
INFILE <- file.dir %&% "results_rasqual_chr" %&% chromo %&% ".txt.gz"
OUTFILE <- file.dir %&% "results_rasqual_chr" %&% chromo %&% "_with_Pvalues.txt"
#INFILE<-paste("results_rasqual_chr",paste(chromo,".txt.gz",sep=""),sep="")
#OUTFILE<-paste("results_rasqual_chr",paste(chromo,"_with_Pvalues.txt",sep=""),sep="")
data <- fread(INFILE,header=FALSE,sep="\t")
names(data)<-c("feature_id","snp_chrpos","chr","snp_pos","ref_allele","alt_allele","allele_freq","hwe_chi2","imputation_score","log10_benjamini_hochberg_qvalue","chi2","pi","delta","phi","overdispersion","snp_number","n_fSNPs","n_rSNPs","n_iterations_H0","n_iterations_Ha","random_location_ties","log_likelihood_H0","convergence_status","R2_fSNPs","R2_rSNPs")
data$p_chi2<-pchisq(data$chi2, df=1, lower.tail=FALSE)
peaks<-unique(data$feature_id)
data$p_bonf_adj<-NA
print("")
print("Calculating bonferroni-adjusted p-values")
pb <- txtProgressBar(min=0,max=length(peaks),style=3)
for (i in 1:length(peaks)){
	setTxtProgressBar(pb,i)
	PEAK<-peaks[i]
	#print(paste(i,PEAK,sep=" "))
	N_TESTS<-max(data$n_rSNPs[data$feature_id==PEAK])
	#print(N_TESTS)
	data$p_bonf_adj[data$feature_id==PEAK]<-p.adjust(data$p_chi2[data$feature_id==PEAK],method="bonferroni")
}
write.table(data,OUTFILE,sep="\t",quote=F,row.names=F)
print("done")
