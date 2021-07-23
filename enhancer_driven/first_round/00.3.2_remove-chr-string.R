"%&%" <- function(a,b) paste0(a,b)
library("data.table");library("dplyr");library("purrr")
work.dir <- "/well/mccarthy/users/jason/projects/wtsa/enhancer_driven/first_round/"
oligo.df <- fread(work.dir %&% "input_files/Oligo-file_merged_GRCh38.txt",
                  header=F)
oligo.df$V2 <- purrr::map(oligo.df$V2,function(s){gsub("chr","",s)}) %>%
as.character(.)
oligo.df$V5 <- purrr::map(oligo.df$V5,function(s){gsub("chr","",s)}) %>%
as.character(.)
write.table(x=oligo.df,
  file=work.dir%&%"input_files/Oligo-file_merged_GRCh38_no-chr.txt",
  quote=F,sep="\t",row.names=F,col.names=F)
