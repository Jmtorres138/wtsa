"%&%" <- function(a,b) paste0(a,b)
library(rtracklayer)
library("data.table")
library("tidyverse")
rescomp.dir <- "/Users/jasont/science/servers/FUSE5/"
work.dir <- rescomp.dir %&% 
  "projects/wtsa/promoter_driven/"
input.dir <- work.dir %&% "input_files/"
gtf.dir <- rescomp.dir %&% "datasets/GTEx/TOPMED_reference_files/"

gtf.file <- gtf.dir %&% "gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf.gz"
gtf <- rtracklayer::import(gtf.file)
gtf.df <- as.data.frame(gtf)

gtf.file2 <- gtf.dir %&% "gencode.v30.GRCh38.annotation.ERCC.gtf.gz"
gtf2 <- rtracklayer::import(gtf.file2)
gtf.df2 <- as.data.frame(gtf2)

gtf.df %>% filter(type=="transcript") %>% dim(.)
gtf.df2 %>% filter(type=="transcript") %>% dim(.)
trans.df <- gtf.df2 %>% filter(type=="transcript")

# Will write a bed file with all TRANSCRIPTs to use a reference in plots 
bed.df <- data.frame(V1=trans.df$seqnames,V2=trans.df$start,V3=trans.df$end,
                     V4=trans.df$gene_name,stringsAsFactors = F)
write.table(x=bed.df,file=input.dir%&%
      "gencode.v30.GRCh38.annotation.ERCC.transcripts.bed",
      sep="\t",quote=F,row.names=F,col.names=F)


