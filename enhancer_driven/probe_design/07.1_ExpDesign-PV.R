
"%&%" <- function(a,b) paste0(a,b) 
library("data.table")
library("dplyr")
library("GenomicRanges")
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library("BSgenome")

#serv.dir <- "/Users/jtorres/FUSE/"
serv.dir <- "/well/got2d/jason/"
work.dir <- serv.dir %&% "projects/wtsa/enhancer_driven/probe_design/"
rds.dir <- work.dir %&% "rds/"
profile.dir <- work.dir %&% "profile-snps/"

prof.df <- fread(profile.dir%&%"profile_credt2d.txt")

save.dir <- work.dir %&% "experimental-design-pv/"
dir.create(save.dir)

gread <- function(bed.df){
  # names must be c(chr, start, end, id, score, strand), where available 
  if(length(bed.df)==3){
    gr <- with(bed.df, GRanges(chr, IRanges(start, end)))
  } else if (length(bed.df)==4){
    gr <- with(bed.df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(bed.df)==5){
    gr <- with(bed.df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(bed.df)>=6){
    bed.df$strand <- gsub(pattern="[^+-]+",replacement='*', bed.df$strand)
    gr <- with(bed.df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}


cred_extract <- function(segnum){
  temp <- filter(prof.df,SEGNUMBER==segnum) %>% arrange(desc(PPA.fgwas))
  hp <- temp$POS[1] # high-priority variant 
  cc <- arrange(temp,desc(PPA))$POS[1] # credible set control 
  cc <- ifelse(cc!=hp,cc,NA) 
  ac <- filter(temp,!(POS %in% c(hp,cc)),endo.atac==TRUE)$POS[1] # atac control 
  #in credible set 
  sub <- filter(temp,POS %in% c(hp,cc,ac))
  expname <- as.character(sapply(1:dim(sub)[1], function(i){
    pos <- sub$POS[i]
    en <- ifelse(pos==hp,"high.priority",
                 ifelse(pos==cc,"cred.control",
                        ifelse(pos==ac,"atac.control",NA)))
    return(en)
  }))
  out.df <- cbind(sub,expname)
  return(out.df)
}


imr90 <- fread(serv.dir%&%
                 "reference/tads/IMR90/combined/total.combined.domain")
hesc <- fread(serv.dir%&%
                "reference/tads/hESC/combined/total.combined.domain")
cell <- c(rep("IMR90",dim(imr90)[1]),rep("hESC",dim(hesc)[1]))
tad.df <- as.data.frame(cbind(rbind(imr90,hesc),cell))
names(tad.df) <- c("chr","start","end","id")
tad.df$start <- tad.df$start + 1
tad.df$end <- tad.df$end + 1
tad.gr <- gread(tad.df)


atac.df <- readRDS(rds.dir%&%"atac.df.RDS")
atac.df$chr <- gsub("chr","ch",atac.df$chr)
atac.gr <- gread(atac.df)



get_external_control <- function(segnum,window=5000){
  temp <- filter(prof.df,SEGNUMBER==segnum) %>% arrange(desc(PPA.fgwas))
  chrom <- temp$CHR[1]
  mn <- min(temp$POS); mx <- max(temp$POS)
  cred.gr <- GRanges(chrom,IRanges(mn,mx))
  tad.ovlp <- tad.gr[tad.gr %over% cred.gr]
  tad.constrain <- GRanges(chrom,IRanges(min(start(tad.ovlp)),max(end(tad.ovlp))))
  
  mymin <- mn - window; mymax <- mx + window 
  chrom <- gsub("chr","ch",chrom)
  min.gr <- GRanges(chrom,IRanges(start(tad.constrain),mymin))
  max.gr <- GRanges(chrom,IRanges(mymax,end(tad.constrain)))
  upstream.snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37,ranges = min.gr)
  downstream.snps <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP144.GRCh37,ranges = max.gr)
  # Subset to ATAC overlaps 
  up.atac.snps <- upstream.snps[upstream.snps %over% atac.gr]
  down.atac.snps <- downstream.snps[downstream.snps %over% atac.gr]
  up <- max(pos(up.atac.snps)); down <- min(pos(down.atac.snps))
  up.dist <- abs(mymin-up); down.dist <- abs(down-mymax)
  up.rs <- up.atac.snps$RefSNP_id[which(pos(up.atac.snps)==up)]
  down.rs <- down.atac.snps$RefSNP_id[which(pos(down.atac.snps)==down)]
  pos <- ifelse(up.dist<down.dist,up,down)
  rsid <- ifelse(up.dist<down.dist,up.rs,down.rs)
  chrom <- temp$CHR[1]
  return(c(chrom,pos,rsid))
}


build_stack <- function(seg){
  temp <- cred_extract(seg)
  temp2 <- temp[1,]
  extcon <-  get_external_control(seg)
  temp2$SNPID <- extcon[3]
  temp2$POS <- extcon[2] 
  temp2$endo.atac <- TRUE 
  temp2[1,c(5,6,7,10,11)] <- NA
  temp2$expname <- "external.atac.control"
  out.df <- rbind(temp,temp2)
  return(out.df)
}

build_full <- function(){
  out.df <- c()
  pb <- txtProgressBar(min=0,max=length(unique(prof.df$SEGNUMBER)),style=3)
  for (i in 1:length(unique(prof.df$SEGNUMBER))){
    setTxtProgressBar(pb,i)
    segnum <- unique(prof.df$SEGNUMBER)[i]
    print(paste("index: ",i, " segnum: ",segnum))
    sub <- filter(prof.df,SEGNUMBER==segnum)
    chr <- sub$CHR[1]
    locus <- sub$LOCUS[1]
    locus <- gsub("/","-",locus); locus <- gsub(" ","-",locus)
    locus <- gsub("(","",locus,fixed = TRUE)
    locus <- gsub(")","",locus,fixed = TRUE)
    savename <- save.dir %&% segnum %&% "_" %&% chr %&% "_" %&% locus %&% ".txt"
    stack.df <- build_stack(segnum)
    write.table(x=stack.df,file=savename,sep="\t",row.names=F,
                quote=F)
    out.df <- rbind(out.df,stack.df)
  }
  write.table(x=out.df,file=save.dir%&%"design-pv.txt",
              sep="\t",row.names=F,quote=F)
}

build_full() 


