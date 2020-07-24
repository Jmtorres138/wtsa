"%&%" <- function(a,b) paste0(a,b)
library("data.table");library("dplyr")
serv.dir <- "/well/mccarthy/users/jason/" 
serv.dir0 <- "/well/got2d/jason/" 
work.dir <- serv.dir %&% "projects/wtsa/joint_analyses/"
coloc.dir <- work.dir %&% "coloc_analysis_files/"
res.dir <- coloc.dir %&% "gwas-t2d_islet-eqtl_coloc/"
out.dir <- coloc.dir %&% "summary_files/"
ref.dir <- serv.dir %&% "datasets/"
gene.ref.file <- ref.dir %&% "gencode.v19.genes.v7.patched_contigs.gtf"

gtf.df <- import(gene.ref.file) %>% as.data.frame()

collate_coloc_results <- function(){
  coloc.sum.df <- c()
  coloc.sig.df <- c()
  rds.vec <- list.files(res.dir)[1:10] # TESTING FIRST 100
  pb <- txtProgressBar(min=0,max=length(rds.vec),style=3)
  for (i in 1:length(rds.vec)){
    setTxtProgressBar(pb,i)
    rds.file <- rds.vec[i]
    s.vec <- strsplit(x=rds.file,split=":")[[1]]
    gwas.signal <- s.vec[1]
    egene <- strsplit(s.vec[2],split=".coloc.RDS")[[1]]
    rds <- readRDS(res.dir%&%rds.file)
    gene.name <- (filter(gtf.df,gene_id==egene)$gene_name %>% unique(.))[1]
    gene.name <- ifelse(length(gene.name)>0,gene.name,NA)
    name.df <- data.frame("gwas.signal"=gwas.signal,"egene"=egene,
                          "egene.name"=gene.name,stringsAsFactors = FALSE)
    build.df1 <- rds$summary %>% t(.) %>% as.data.frame(.)
    build.df2 <- rds$results %>% arrange(.,desc(SNP.PP.H4)) %>% 
      filter(.,SNP.PP.H4>=0.5)     
    build.df1 <- cbind(name.df,build.df1)
    coloc.sum.df <- rbind(coloc.sum.df,build.df1)
    if (dim(build.df2)[1]>0){
      build.df2 <- cbind(name.df,build.df2)
      coloc.sig.df <- rbind(coloc.sig.df,build.df2)
    }
  }
  return(list(coloc.sum.df,coloc.sig.df))
}

coloc.list <- collate_coloc_results()
write.table(x=coloc.list[[1]],file=out.dir %&%
              "coloc_t2d-gwas-islet-eqtl_summary.txt",
            sep="\t",quote=F,row.names=F)
write.table(x=coloc.list[[2]],file=out.dir %&%
              "coloc_t2d-gwas-islet-eqtl_sig-p50.txt",
            sep="\t",quote=F,row.names=F)
