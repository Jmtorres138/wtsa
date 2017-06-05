
"%&%" <- function(a,b) paste0(a,b) 
library("data.table")
library("dplyr")
library("GenomicRanges")
library(Homo.sapiens)
library("ggbio")
library("ggrepel")
serv.dir <- "/Users/jtorres/FUSE/"
work.dir <- serv.dir %&% "projects/wtsa/enhancer_driven/fgwas_probe_design/first_round/"
rds.dir <- work.dir %&% "RDS/"
ld.dir <- work.dir %&% "ld/ld_files/"
all.eqtl.dir <- serv.dir %&% "reference/islet/eqtls/inspire/nominal_pass/output/"


credt2d.df <- readRDS(rds.dir%&%"credComb.1000less.df.RDS")
credt2d.df$CHR <- "chr" %&% credt2d.df$CHR
eqtl.df <- readRDS(rds.dir%&%"eqtl.fdr01.df.RDS")
atac.df <- readRDS(rds.dir%&%"atac.df.RDS")
gcred.df <- readRDS(rds.dir%&%"genetic-credible-sets-ind.RDS")



dpn.file <- serv.dir %&% "reference/DpnII/hg19_DpnII-sites.bed.gz" #  0-based bed file 
dpn.df <- fread("cat " %&% dpn.file %&% " | zmore")
# adjust to 1-based scheme 
dpn.df$V2 <- dpn.df$V2 + 1 
dpn.df$V3 <- dpn.df$V3 + 1 



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

get_overlaps <- function(a,b){
  # a and b are genomic ranges
  vec <- a %over% b
  ##print(sum(vec) %&% " / " %&% length(vec))
  return(vec)
}



# Plotting Functions


fcredld_plot <- function(segnum,fac=10){
  df <- filter(credt2d.df,SEGNUMBER==segnum) %>% arrange(desc(PPA))
  mx <- max(df$POS); mn <- min(df$POS); span <- mx - mn; loc <- df$refGene[1]
  span <- ifelse(span==0,1,span)
  chrom <- df$CHR[1];mymin <- (mn-fac*span);mymax <- (mx+fac*span)  
  loc <- gsub("/","-",loc)
  expan <- 1e6 
  f <- chrom%&%"_"%&%(mn-expan)%&%"-"%&%(mx+expan)%&%"_"%&%loc%&%"_"%&%segnum%&%".hap.ld"
  temp.df <- fread(ld.dir %&% f)
  temp.df <- temp.df[,1:7]
  names(temp.df) <- c("CHR","POS1",	"POS2","N_CHR","R2","D","Dprime")
  maxpos <- df$POS[1] # by PPA 
  LeadR2 <- as.numeric(sapply(1:length(df$POS),function(i){
    pos <- df$POS[i]
    if (pos==maxpos){r2<-1}
    else {r2<-filter(temp.df,(POS1==maxpos & POS2==pos) | 
                       (POS1==pos & POS2==maxpos))$R2}
    return(r2)
  }))
  Color <- as.character(sapply(1:length(LeadR2),function(i){
    r2 <- LeadR2[i]
    col <- ifelse(r2>=0.80,"red", 
                  ifelse(r2<0.8 & r2>=0.6, "yellow",
                         ifelse(r2<0.6 & r2>=0.4, "green",
                                ifelse(r2<0.4 & r2>=0.2, "dodgerblue1",
                                       ifelse(r2<0.2 & r2>=0.0, "purple4",
                                              "purple4")))))
  }))
  df <- cbind(df,LeadR2,Color)
  plt <- ggplot(data=df,aes(x=POS,y=PPA)) + 
    geom_point(shape=21,color="black",size=2,fill=df$Color) + 
    scale_fill_continuous(low="yellow",
                          high="firebrick1",limits=c(0,1),
                          breaks=seq(0,1,0.2)) + 
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1)) + 
    theme_bw() + 
    theme(legend.position="none",
          panel.grid.minor=element_blank(),
          axis.title.y=element_text(size=5),
          axis.text.y=element_text(size=4)) +
    scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1)) + 
    coord_cartesian(ylim=c(0,1),expand=TRUE)
  plt <- plt + geom_label_repel(data=arrange(df,desc(PPA))[1:1,],
                                aes(x=POS,y=PPA,label=SNPID),size=3)
  return(plt)
}


gcred_plot <- function(segnum,fac=10){
  df <- filter(credt2d.df,SEGNUMBER==segnum) %>% arrange(desc(PPA))
  mx <- max(df$POS); mn <- min(df$POS); span <- mx - mn; loc <- df$refGene[1]
  span <- ifelse(span==0,1,span)
  chrom <- df$CHR[1];mymin <- (mn-fac*span);mymax <- (mx+fac*span)  
  
  gloc <- filter(gcred.df,CHR==(gsub("chr","",chrom)),POS>=mn,POS<=mx)$LOCUS
  sub.df <- filter(gcred.df,LOCUS==gloc)
  plt <- ggplot(data=sub.df,aes(x=POS,y=PPA)) + 
    geom_point(shape=21,color="black",size=2,fill="grey") + 
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1)) + 
    theme_bw() + 
    theme(legend.position="none",
          panel.grid.minor=element_blank(),
          axis.title.y=element_text(size=5),
          axis.text.y=element_text(size=4)) +
    scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1)) + 
    coord_cartesian(ylim=c(0,1),expand=TRUE)
  plt <- plt + geom_label_repel(data=arrange(sub.df,desc(PPA))[1:1,],
                                aes(x=POS,y=PPA,label=SNPID),size=3)
  return(plt)
}


atac_plot <- function(segnum,fac=0.10){
  df <- filter(credt2d.df,SEGNUMBER==segnum)
  mx <- max(df$POS); mn <- min(df$POS); span <- mx - mn
  span <- ifelse(span==0,1,span)
  chrom <- df$CHR[1];mymin <- (mn-fac*span);mymax <- (mx+fac*span)
  a1sub <- filter(atac.df,chr==chrom,start>=mymin,end<=mymax,id=="Oxford") %>% 
    dplyr::select(one_of("chr","start","end"))
  a1.g <- gread(a1sub[,1:3])
  a2sub <- filter(atac.df,chr==chrom,start>=mymin,end<=mymax,id=="Parker") %>% 
    dplyr::select(one_of("chr","start","end"))
  a2.g <- gread(a2sub[,1:3])
  if (length(a1.g)>0 & length(a2.g)>0){
    plt <- ggplot(a2.g) + geom_segment(color="firebrick3",size=30,alpha=0.9) +
      geom_segment(a1.g,color="yellow2",size=25,alpha=0.9) + ylab("") +
      theme_bw() + theme(panel.grid = element_blank())
  } else if (length(a1.g)>0){
    plt <- ggplot(a1.g) + geom_segment(color="yellow2",size=25,alpha=0.9) +
      theme_bw() + theme(panel.grid = element_blank())
  } else if (length(a2.g)>0){
    plt <- ggplot(a2.g) + geom_segment(color="firebrick3",size=30,alpha=0.9) +
      theme_bw() + theme(panel.grid = element_blank())
  } else {
    dframe <- data.frame(start=mymin,end=mymax)
    plt <- ggplot(data=dframe) + geom_rect(data=dframe,aes(xmin=start,xmax=end),
                                           ymin=0,ymax=0) + theme_bw()
  }
  return(plt)
}



dpnII_plot <- function(segnum,fac=0.10){
  df <- filter(credt2d.df,SEGNUMBER==segnum)
  mx <- max(df$POS); mn <- min(df$POS); span <- mx - mn
  span <- ifelse(span==0,1,span)
  chrom <- df$CHR[1];mymin <- (mn-fac*span);mymax <- (mx+fac*span)
  dpnsub <- filter(dpn.df,V1==chrom,V2>=mymin,V3<=mymax)
  if (dim(dpnsub)[1]==0){
    dpnsub <- data.frame(V1=chrom,V2=mymin,V3=mymax)
    plt <-  ggplot(data=dpnsub)
    plt <- plt + geom_vline(aes(xintercept=V2), size=0,color="white") + 
      geom_vline(aes(xintercept=V3),size=0,color="white")
  } else{
    plt <-  ggplot(data=dpnsub)
    plt <- plt + geom_vline(aes(xintercept=V2), size=0.1) + 
      geom_vline(aes(xintercept=V3),size=0.1) 
  }
  plt <- plt + ylim(c(0,1)) + 
    theme_bw() + theme(axis.text.y=element_blank(),
                       axis.ticks = element_blank(),
                       panel.grid = element_blank())  
  
  return(plt)
}



gene_plot <- function(segnum,fac=0.10){
  df <- filter(credt2d.df,SEGNUMBER==segnum)
  mx <- max(df$POS); mn <- min(df$POS)
  chrom <- df$CHR[1]; span <- mx - mn
  span <- ifelse(span==0,1,span)
  gr <- GRanges(chrom,IRanges(start=mn-1e6,end=mx+1e6))
  aplt <- autoplot(Homo.sapiens, which = gr, label.color = "black",#
                   color = "brown",
                   #fill = "brown") + #,stat="reduce"
                   fill = "brown") + #,stat="reduce") + 
    xlim(c(mn,mx)) + 
    theme_alignment(grid = FALSE,border = FALSE) 
  return(aplt)
}


islet_eqtl_plot <- function(segnum,fac=0.10){
  df <- filter(credt2d.df,SEGNUMBER==segnum)
  mx <- max(df$POS); mn <- min(df$POS)
  chrom <- df$CHR[1]; loc <- df$LOCUS[1]; span <- mx - mn
  span <- ifelse(span==0,1,span)
  mymin <- (mn-fac*span); mymax <- (mx+fac*span)
  sub.df <- filter(eqtl.df, CHR==chrom, POS>=mymin) %>% filter(POS<=mymax)
  maxy <- max(-log(sub.df$P,base=10))
  maxy <- ifelse(maxy<10,10,maxy+2)
  plt <- ggplot(data=sub.df, aes(x=POS,y=-log(P,base=10))) +  
    geom_point(shape=21,color="black", size=2,
               fill="green1",alpha=0.80) + 
    ylab(expression(paste("-log"[10],"(p-value)"))) + theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.title.y=element_text(size=5),
          axis.text.y=element_text(size=4)) +  xlim(c(mn,mx)) +
    ylim(c(0,maxy))
  return(plt)
}


prob.dir <- serv.dir %&% "projects/wtsa/enhancer_driven/" %&% "probe_design/from_Martijn/DpnII_120bp/"
probe.df <- fread(prob.dir %&% "lessThan100_all_success.txt")
chrom <- as.character(sapply(1:dim(probe.df)[1],function(i){
  frag <- probe.df$`Fragment Coordinates`[i]
  c <- strsplit(x=frag,split=":")[[1]][1]
  return(c)
}))
frag.start <- as.integer(sapply(1:dim(probe.df)[1],function(i){
  frag <- probe.df$`Fragment Coordinates`[i]
  span <- strsplit(x=frag,split=":")[[1]][2]
  fs <- as.integer(strsplit(x=span,split="-")[[1]][1])
  return(fs)
}))
frag.end <- as.integer(sapply(1:dim(probe.df)[1],function(i){
  frag <- probe.df$`Fragment Coordinates`[i]
  span <- strsplit(x=frag,split=":")[[1]][2]
  fe <- as.integer(strsplit(x=span,split="-")[[1]][2])
  return(fe)
}))
left.start <- as.integer(sapply(1:dim(probe.df)[1],function(i){
  span <- probe.df$`Left Oligo`[i]
  ls <- as.integer(strsplit(x=span,split="-")[[1]][1])
  return(ls)
}))
left.end <- as.integer(sapply(1:dim(probe.df)[1],function(i){
  span <- probe.df$`Left Oligo`[i]
  le <- as.integer(strsplit(x=span,split="-")[[1]][2])
  return(le)
}))
right.start <- as.integer(sapply(1:dim(probe.df)[1],function(i){
  span <- probe.df$`Right Oligo`[i]
  rs <- as.integer(strsplit(x=span,split="-")[[1]][1])
  return(rs)
}))
right.end <- as.integer(sapply(1:dim(probe.df)[1],function(i){
  span <- probe.df$`Right Oligo`[i]
  re <- as.integer(strsplit(x=span,split="-")[[1]][2])
  return(re)
}))
probe.df <- cbind(probe.df,chrom,frag.start,frag.end,
                  left.start,left.end,right.start,right.end)

prob_plot <- function(segnum,fac=0.1){
  df <- filter(credt2d.df,SEGNUMBER==segnum)
  mx <- max(df$POS); mn <- min(df$POS)
  mychrom <- df$CHR[1]; loc <- df$LOCUS[1]; span <- mx - mn
  span <- ifelse(span==0,1,span)
  mymin <- (mn-fac*span); mymax <- (mx+fac*span)
  sub <- filter(probe.df,chrom==mychrom)
  if (dim(sub)[1]==0){
    sub <- data.frame(a=NA,b=NA,c=NA,d=NA,e=NA,f=NA,g=NA,h=NA,i=NA,j=NA,k=NA,l=NA)
    names(sub) <- names(probe.df)
  }
  plt <- ggplot(data=sub) + 
    geom_rect(sub, fill="black", aes(xmin=frag.start,xmax=frag.end),
              ymin=0.4, ymax=0.6) + 
    geom_rect(sub, fill="dodgerblue2", aes(xmin=left.start,xmax=left.end),
              ymin=0.6, ymax=0.7) +
    geom_rect(sub, fill="dodgerblue2", aes(xmin=right.start,xmax=right.end),
              ymin=0.3, ymax=0.4) +
    coord_cartesian(xlim=c(mymin,mymax)) + theme_bw() + 
    theme(panel.grid = element_blank())
  return(plt)
}





# Track Plot function 



track_plot <- function(segnum){
  df <- filter(credt2d.df,SEGNUMBER==segnum)
  mx <- max(df$POS); mn <- min(df$POS); span <- mx - mn
  span <- ifelse(span==0,1,span)
  fac <- ifelse(span==1,10000, 
                ifelse(span>1&span<10000,1,
                       ifelse(span>10000,0.1,0.1)))
  mymin <- (mn-fac*span); mymax <- (mx+fac*span)
  loc <- df$refGene[1]
  gloc <- filter(gcred.df,CHR==(gsub("chr","",chrom)),POS>=mn,POS<=mx)$LOCUS
  loc <- ifelse(loc==gloc,loc,(loc%&%"/"%&%gloc))
  atac_plt <- atac_plot(segnum,fac)
  dpn_plt <- dpnII_plot(segnum,fac)
  fcredld_plt <- fcredld_plot(segnum,fac)
  gene_plt <- gene_plot(segnum,fac)
  probe_plt <- prob_plot(segnum,fac)
  eqtl_plt <- islet_eqtl_plot(segnum,fac)
  gcred_plt <- gcred_plot(segnum,fac)
  
  tracks(`Islet eQTL \n(FDR<0.1)`=eqtl_plt,
         `ATAC \nEndoC`=atac_plt,
         `Functional Credible Set \n(LD)`=fcredld_plt,
         `Genetic Credible Set`=gcred_plt,
         `Oligos`=probe_plt,
         `DpnII`=dpn_plt,
         `Genes`=gene_plt,
         heights=c(0.5,0.5,2,2,0.2,0.5,2),title=loc %&% "\nInterval: " %&% 
           
           round(span/1000,digits=2) %&%" kb",
         label.text.cex=0.6,main.height=2) +
    scale_x_sequnit("Mb") + 
    theme(axis.text.x=element_text(size=5)) 
}





