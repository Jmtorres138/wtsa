"%&%" <- function(a,b) paste0(a,b)
library("data.table");library("dplyr")
library("GenomicRanges");library("purrr")
library("coloc")
serv.dir <- "/well/mccarthy/users/jason/" # "/home/jason/science/servers/FUSE5/" # 
serv.dir0 <- "/well/got2d/jason/" # "/home/jason/science/servers/FUSE/" # 
work.dir <- serv.dir %&% "projects/wtsa/joint_analyses/"
coloc.dir <- work.dir %&% "coloc_analysis_files/"
gwas.stat.dir <- coloc.dir %&% "conditioned_gwas_summarystats/"
out.dir <- coloc.dir %&% "gwas-t2d_islet-eqtl_coloc/"
eqtl.dir <- serv.dir0 %&% "reference/islet/eqtls/oxford/nominal_pass/output/"
eqtl.chrom.dir <- eqtl.dir %&% "eqtls_all_by_chrom/"

args <- commandArgs(trailingOnly=TRUE)
chromo <- args[1]

eqtl.df <- fread(eqtl.chrom.dir %&% "eqtl_all_" %&% chromo %&% ".txt.gz")
names(eqtl.df) <- c("gene.id","snp.id","dist","nom.p","nom.beta","chrom",
                    "pos","ref","alt")
eqtl.df$pos <- as.integer(eqtl.df$pos)
eqtl.df$chrom <- "chr" %&% eqtl.df$chrom
eqtl.gr <- GRanges(seqnames=unique(eqtl.df$chrom),
                   IRanges(start=eqtl.df$pos,end=eqtl.df$pos))


file.vec <- list.files(gwas.stat.dir)
keep.vec <- purrr::map(file.vec,function(s){
  grepl(pattern="-"%&%chromo%&%"-",x=s)
}) %>% as.logical(.)
file.vec <- file.vec[keep.vec]

gwas_eqtl_coloc <- function(gwas.signal.dir,gwas.signal,
                            n.gwas=898130,case.prop=0.09, n.eqtl=174){
  gwas.signal.file <- gwas.signal.dir %&% gwas.signal %&% ".csv" 
  # Extract all eQTL info for all SNPs that overlap GWAS signal 
  sig.df <- fread(gwas.signal.file)
  s <- strsplit(gwas.signal.file,split=".",fixed=TRUE)[[1]][[1]]
  vec <- strsplit(s,split="-")[[1]]
  pos <- vec[length(vec)] %>% as.integer(.)
  region.gr <- GRanges(seqnames = unique(sig.df$CHR),
                       IRanges(start=(pos-500000),end=(pos+500000)))
  eqtl.sub.df <- eqtl.df[(eqtl.gr %over% region.gr),]
  gene.vec <- eqtl.sub.df$gene.id %>% unique(.)
  pb <- txtProgressBar(min=1,max=length(gene.vec),style=3)
  for (i in 1:length(gene.vec)){
    #print(i)
    setTxtProgressBar(pb,i)
    gene <- gene.vec[i]
    gene.sub.df <- filter(eqtl.sub.df,gene.id==gene)
    sig.sub.df <- sig.df[(sig.df$POS %in% gene.sub.df$pos),]
    overlap.pos <- sig.sub.df$POS
    overlap.snps <- sig.sub.df$SNPID
    names(overlap.pos) <- overlap.snps
    p.vec.gwas <- c()
    p.vec.eqtl <- c() 
    maf.vec <- c()
    for (posit in overlap.pos){
      p1 <- filter(sig.df,POS==posit)$PVAL
      names(p1) <- filter(sig.df,POS==posit)$SNPID
      p2 <- filter(gene.sub.df,pos==posit)$nom.p
      names(p2) <- filter(sig.df,POS==posit)$SNPID
      f <- filter(sig.df,POS==posit)$'F'
      maf <- min(f,(1-f))
      p.vec.gwas <- append(p.vec.gwas,p1)
      p.vec.eqtl <- append(p.vec.eqtl,p2)
      maf.vec <- append(maf.vec,maf)
    }
    exclude.vec <- is.na(maf.vec) | is.na(p.vec.gwas) | is.na(p.vec.eqtl) 
    p.vec.gwas <- p.vec.gwas[!exclude.vec]
    p.vec.eqtl <- p.vec.eqtl[!exclude.vec]
    maf.vec <- maf.vec[!exclude.vec]
    if (length(maf.vec)>0 & length(maf.vec)==length(p.vec.gwas) & 
        length(maf.vec==p.vec.eqtl)){
        coloc.res <- suppressMessages(coloc.abf(dataset1=list(pvalues=p.vec.gwas,
                                             snp=names(p.vec.gwas),
                                           N=n.gwas,s=case.prop,type="cc"),
                             dataset2=list(pvalues=p.vec.gwas,N=n.eqtl,
                                           snp=names(p.vec.gwas),type="quant"),
                             MAF=maf.vec))
        out.name <- out.dir %&% gwas.signal %&% ":" %&% gene %&% ".coloc.RDS"
        saveRDS(object=coloc.res,file=out.name)      
    } else{
      print("Not enough information present for Gene: " %&% gene)
    }
  }
}

for (i in 1:length(file.vec)){
  gwas.signal <- strsplit(file.vec[i],split=".",fixed=TRUE)[[1]][1]
  print(gwas.signal)
  gwas_eqtl_coloc(gwas.stat.dir,gwas.signal)  
}



