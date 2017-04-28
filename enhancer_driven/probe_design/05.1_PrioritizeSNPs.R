#For each SNP at each locus will determine if it (1) overlaps an EndoC-BH1 ATAC region, (2) PPA improves with islet annotations i.e. fgwas NOTE: this is the `change` entry in credt2d.df, (3) is in LD with islet eQTL


"%&%" <- function(a,b) paste0(a,b)
library("data.table")
library("dplyr")

serv.dir <- "/well/got2d/jason/"
rds.dir <- serv.dir %&% "projects/wtsa/enhancer_driven/probe_design/rds/"
source(serv.dir %&% "projects/wtsa/enhancer_driven/probe_design/" %&%
         "04.1b_VizTracks-T2D.R")

ld.dir <- serv.dir %&% "projects/wtsa/enhancer_driven/probe_design/ld/ld_files/"
txt.dir <- serv.dir %&% "projects/wtsa/enhancer_driven/probe_design/profile-snps/"

print(rds.dir)


# EVALUATE EndoC-BH1 ATAC peak


eval_endo_atac <- function(segnum,cred.df=credt2d.df){
  sub.df <- filter(cred.df,SEGNUMBER==segnum)
  chrom <- sub.df$CHR[1]
  endo.atac  <- as.logical(sapply(1:dim(sub.df)[1], function(i){
    pos <- sub.df$POS[i];temp <- filter(atac.df,chr==chrom,start<=pos,end>=pos)
    return(ifelse(dim(temp)[1]>0,TRUE,FALSE))
  }))
  sub.df <- cbind(sub.df,endo.atac)
  return(sub.df)
}




#Evalute LD with islet eQTL

eval_islet_eqtl <- function(segnum,cred.df=credt2d.df,ld.thresh=0.8){
  sub.df <- filter(cred.df,SEGNUMBER==segnum)
  chrom <- sub.df$CHR[1]
  mn <- min(sub.df$POS); mx <- max(sub.df$POS); 
  chrom<-sub.df$CHR[1]; locus <- sub.df$LOCUS[1]; study <- sub.df$STUDY[1]
  locus <- gsub("/","-",locus)
  locus <- gsub(" ","-",locus)
  locus <- gsub("(","",locus,fixed=TRUE); locus <- gsub(")","",locus,fixed=TRUE)
  window <- 1e6
  mymin <- mn-window; mymax<-mx+window
  
  eqtl.sub <- filter(eqtl.df,CHR==chrom,POS>=mymin,POS<=mymax)
  
  fname <- chrom%&%"_"%&%mymin%&%"-"%&%mymax%&%"_"%&%locus%&%"_"%&%study%&%".hap.ld"  
  ld.df <- fread(ld.dir%&%fname,select=c(1:3,5))
  names(ld.df) <- c("chr","pos1","pos2","r2")
  
  num.eqtls  <- as.integer(sapply(1:dim(sub.df)[1], function(i){
    pos <- sub.df$POS[i]
    temp <- filter(ld.df,chr==gsub("chr","",chrom)) %>% 
      filter(pos1==pos | pos2==pos) %>% filter(r2 >= ld.thresh)
    pos.vec <- unique(c(temp$pos1,temp$pos2))
    eval <- sum(pos.vec %in% eqtl.sub$POS)
    num.eqtl <- ifelse(dim(temp)[1]>0,ifelse(eval>0,eval,0),0)
    return(num.eqtl)
  }))
  egenes  <- as.character(sapply(1:dim(sub.df)[1], function(i){
    pos <- sub.df$POS[i]
    temp <- filter(ld.df,chr==gsub("chr","",chrom)) %>% 
      filter(pos1==pos | pos2==pos) %>% filter(r2 >= ld.thresh)
    pos.vec <- unique(c(temp$pos1,temp$pos2))
    eval <- sum(pos.vec %in% eqtl.sub$POS)
    egene <- ifelse(eval>0, 
                    paste(unique(filter(eqtl.sub,POS %in% pos.vec)$GENE),collapse=","),NA)
    return(egene)
  }))
  out.df <- cbind(num.eqtls,egenes)
  return(out.df) 
}



#Build Prioritization Data Frame

build_df <- function(cred.df=credt2d.df){
  out.df <- c() 
  pb <- txtProgressBar(min=0,max=length(unique(cred.df$SEGNUMBER)),style=3)
  for (i in 1:length(unique(cred.df$SEGNUMBER))){
    setTxtProgressBar(pb,i)
    segnum <- unique(cred.df$SEGNUMBER)[i]
    print(segnum)
    stack.df <- eval_endo_atac(segnum)
    stack.df <- cbind(stack.df,eval_islet_eqtl(segnum,ld.thresh = 0.80))
    names(stack.df)[dim(stack.df)[2]-1] <- names(stack.df)[dim(stack.df)[2]-1] %&% ".ld80"
    names(stack.df)[dim(stack.df)[2]] <- names(stack.df)[dim(stack.df)[2]] %&% ".ld80"
    stack.df <- cbind(stack.df,eval_islet_eqtl(segnum,ld.thresh = 0.20))
    names(stack.df)[dim(stack.df)[2]-1] <- names(stack.df)[dim(stack.df)[2]-1] %&% ".ld20"
    names(stack.df)[dim(stack.df)[2]] <- names(stack.df)[dim(stack.df)[2]] %&% ".ld20"
    write.table(x=stack.df,file=txt.dir%&%"profile_seg"%&%segnum%&%".txt",
                quote = F,sep="\t",row.names=F)
    out.df <- rbind(out.df,stack.df)
  }
  out.df$num.eqtls.ld80 <- as.integer(as.character(out.df$num.eqtls.ld80))
  out.df$egenes.ld80 <- as.character(out.df$egenes.ld80)
  out.df$num.eqtls.ld20 <- as.integer(as.character(out.df$num.eqtls.ld20))
  out.df$egenes.ld20 <- as.character(out.df$egenes.ld20)
  return(out.df)
}

ptzd.df <- build_df()
write.table(x=ptzd.df,file=txt.dir%&%"profile_credt2d.txt",
            quote = F,sep="\t",row.names=F)
