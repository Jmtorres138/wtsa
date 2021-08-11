"%&%" <- function(a,b) paste0(a,b)
library("dplyr");library("data.table")
work.dir <- "/well/mccarthy/users/jason/projects/wtsa/enhancer_driven/second_round/"
oligo.file <- work.dir %&% "input_files/Oligo-file_merged_GRCh38_no-chr.txt"
out.dir <- work.dir %&% "03_CaptureCompare/"
oligo.df <- fread(oligo.file,header=F)

out.df <- data.frame("V1"=oligo.df$V1,"V2"="chr"%&%oligo.df$V2,"V3"=oligo.df$V3,
            "V4"=oligo.df$V4,"V5"=oligo.df$V6,"V6"=oligo.df$V7,
            "V7"=oligo.df$V3-200000,"V8"=oligo.df$V4+200000,"V9"=250,"V10"=5000,
          stringsAsFactors=FALSE)
out.df$V7 <- as.integer(out.df$V7)
out.df$V8 <- as.integer(out.df$V8)
write.table(x=out.df,file=out.dir%&%"capture-compare-input_hg38.txt",
  quote=F,sep="\t",row.names=F,col.names=F)
