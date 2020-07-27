"%&%" <- function(a,b) paste0(a,b)
library("data.table");library("dplyr")
library("GenomicRanges");library("purrr")
library("coloc")
serv.dir <- "/well/mccarthy/users/jason/" # "/home/jason/science/servers/FUSE5/" # 
serv.dir0 <- "/well/got2d/jason/" # "/home/jason/science/servers/FUSE/" # 
work.dir <- serv.dir %&% "projects/wtsa/joint_analyses/"
coloc.dir <- work.dir %&% "coloc_analysis_files/"
gwas.stat.dir <- coloc.dir %&% "conditioned_gwas_summarystats/"
out.dir <- coloc.dir %&% "gwas-t2d_islet-caqtl_coloc/"
caqtl.chrom.dir <- coloc.dir %&% "rasqual_output/"
args <- commandArgs(trailingOnly=TRUE)
chromo <- args[1]

caqtl.df <- fread(caqtl.chrom.dir %&% "results_rasqual_" %&% 
                    chromo %&% ".txt.gz")
							
caqtl.df$snp_pos <- as.integer(caqtl.df$snp_pos)
caqtl.gr <- GRanges(seqnames=unique(caqtl.df$chr),
                   IRanges(start=caqtl.df$snp_pos,end=caqtl.df$snp_pos))


file.vec <- list.files(gwas.stat.dir)
keep.vec <- purrr::map(file.vec,function(s){
  grepl(pattern="-"%&%chromo%&%"-",x=s)
}) %>% as.logical(.)
file.vec <- file.vec[keep.vec]

gwas_caqtl_coloc <- function(gwas.signal.dir,gwas.signal,
                            n.gwas=898130,case.prop=0.09, n.caqtl=177){
  gwas.signal.file <- gwas.signal.dir %&% gwas.signal %&% ".csv" 
  # Extract all caQTL info for all SNPs that overlap GWAS signal 
  sig.df <- fread(gwas.signal.file)
  s <- strsplit(gwas.signal.file,split=".",fixed=TRUE)[[1]][[1]]
  vec <- strsplit(s,split="-")[[1]]
  pos <- vec[length(vec)] %>% as.integer(.)
  region.gr <- GRanges(seqnames = unique(sig.df$CHR),
                       IRanges(start=(pos-500000),end=(pos+500000)))
  caqtl.sub.df <- caqtl.df[(caqtl.gr %over% region.gr),]
  peak.vec <- caqtl.sub.df$feature_id %>% unique(.)
  pb <- txtProgressBar(min=1,max=length(peak.vec),style=3)
  for (i in 1:length(peak.vec)){
    #print(i)
    setTxtProgressBar(pb,i)
    peak <- peak.vec[i]
    peak.sub.df <- filter(caqtl.sub.df,feature_id==peak)
    sig.sub.df <- sig.df[(sig.df$POS %in% peak.sub.df$snp_pos),]
    overlap.pos <- sig.sub.df$POS
    overlap.snps <- sig.sub.df$SNPID
    names(overlap.pos) <- overlap.snps
    p.vec.gwas <- c()
    p.vec.caqtl <- c() 
    maf.vec <- c()
    for (posit in overlap.pos){
      p1 <- filter(sig.df,POS==posit)$PVAL
      names(p1) <- filter(sig.df,POS==posit)$SNPID
      p2 <- filter(peak.sub.df,snp_pos==posit)$p_chi2
      names(p2) <- filter(sig.df,POS==posit)$SNPID
      f <- filter(sig.df,POS==posit)$'F'
      maf <- min(f,(1-f))
      if (length(p1)==length(p2) & length(p2)==length(maf)){
        p.vec.gwas <- append(p.vec.gwas,p1)
        p.vec.caqtl <- append(p.vec.caqtl,p2)
        maf.vec <- append(maf.vec,maf)       
      } else{
        print(posit)
      }
    }
    exclude.vec <- is.na(maf.vec) | is.na(p.vec.gwas) | is.na(p.vec.eqtl) 
    p.vec.gwas <- p.vec.gwas[!exclude.vec]
    p.vec.caqtl <- p.vec.caqtl[!exclude.vec]
    maf.vec <- maf.vec[!exclude.vec]
    if (length(maf.vec)>0 & length(maf.vec)==length(p.vec.gwas) & 
        length(maf.vec==p.vec.eqtl)){
        coloc.res <- suppressMessages(coloc.abf(dataset1=list(pvalues=p.vec.gwas,
                                             snp=names(p.vec.gwas),
                                           N=n.gwas,s=case.prop,type="cc"),
                             dataset2=list(pvalues=p.vec.caqtl,N=n.caqtl,
                                           snp=names(p.vec.caqtl),type="quant"),
                             MAF=maf.vec))
        out.name <- out.dir %&% gwas.signal %&% ":" %&% peak %&% ".coloc.RDS"
        saveRDS(object=coloc.res,file=out.name)      
    } else{
      print("Not enough information present for Peak: " %&% peak)
    }
  }
}

for (i in 1:length(file.vec)){
  gwas.signal <- strsplit(file.vec[i],split=".",fixed=TRUE)[[1]][1]
  print(gwas.signal)
  gwas_caqtl_coloc(gwas.stat.dir,gwas.signal)  
}



