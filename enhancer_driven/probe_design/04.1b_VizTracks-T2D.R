
"%&%" <- function(a,b) paste0(a,b) 
library("data.table")
library("dplyr")
library("GenomicRanges")

serv.dir <- "/Users/jtorres/FUSE/"
rds.dir <- serv.dir %&% "projects/wtsa/enhancer_driven/probe_design/rds/"



credt2d.df <- readRDS(rds.dir%&%"credt2d.df.RDS")
eqtl.df <- readRDS(rds.dir%&%"eqtl.df.RDS")
atac.df <- readRDS(rds.dir%&%"atac.df.RDS")


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


library(ggbio)
library(Homo.sapiens)

viz_loc <- function(segnum,fac=0.10,genes=TRUE){
  df <- filter(credt2d.df,SEGNUMBER==segnum)
  mx <- max(df$POS); mn <- min(df$POS)
  (span <- mx - mn)
  print(span)
  chrom <- df$CHR[1] 
  mymin <- (mn-fac*span)
  mymax <- (mx+fac*span)
  loc <- df$LOCUS[1]
  
  # CREATE DPNII / Credible Set Plot 
  dpnsub <- filter(dpn.df,V1==chrom,V2>mymin,V3<mymax)
  # Credible Set Plot 
  plt1 <- ggplot(df,aes(x=POS)) + 
                geom_point(aes(y=PPA),shape=21,color="gray",fill="gray") + 
                geom_point(aes(y=PPA.fgwas,fill=(change>0)),shape=23,color="black") +   
                scale_fill_manual(values=c("dodgerblue2","firebrick2")) + 
    xlab("Position on Chromosome " %&% chrom) + 
    theme_bw() + scale_y_continuous(breaks=seq(0,1,0.05)) + 
    theme(legend.position="none",
          panel.grid.minor=element_blank())
  # DpnII plot
  plt.dpnII <-  ggplot(data=dpnsub) +
    geom_vline(aes(xintercept=V2), size=0.1) + 
    geom_vline(aes(xintercept=V3),size=0.1) + ylim(c(0,1)) + 
    theme_bw() + theme(axis.text.y=element_blank(),
                       axis.ticks = element_blank(),
                       panel.grid = element_blank())
  # Make Reference Gene Track 
  #gr <- GRanges(chrom,IRanges(start=mymin,end=mymax))
  gr <- GRanges(chrom,IRanges(start=mn-1e6,end=mx+1e6))
  if (genes==TRUE){
    ap <- autoplot(Homo.sapiens, which = gr, label.color = "black",#
                 color = "brown",
                 fill = "brown") + #,stat="reduce"
      xlim(c(mn,mx)) + 
      theme_alignment(grid = FALSE,border = FALSE) 
      
  }
  # Make Significant eQTL plot 
  e.sub <- filter(eqtl.df, CHR==chrom, POS >=mymin) %>% filter(POS<=mymax)
  plt2 <- ggplot(data=e.sub, aes(x=POS,y=-log(P,base=10))) +
                   geom_point(color="black",shape=25,fill="green1",size=2) +
    ylab(expression(paste("-log"[10],"(p-value)"))) + theme_bw() + 
    theme(panel.grid = element_blank())
  

  a1sub <- filter(atac.df,chr==chrom,start>=mymin,end<=mymax,id=="Oxford") %>% 
    dplyr::select(one_of("chr","start","end"))
  a1.g <- gread(a1sub[,1:3])
  a2sub <- filter(atac.df,chr==chrom,start>=mymin,end<=mymax,id=="Parker") %>% 
    dplyr::select(one_of("chr","start","end"))
  a2.g <- gread(a2sub[,1:3])
  plt3 <- ggplot(a2.g) + geom_segment(color="firebrick3",size=30,alpha=0.9) +
    geom_segment(a1.g,color="yellow2",size=25,alpha=0.9) + ylab("") +
    theme_bw() + theme(panel.grid = element_blank())
  if (genes==TRUE){
    if(dim(e.sub)[1]>0){
      tracks(`ATAC \nEndoC`=plt3,
             `Islet eQTL\n(FDR < 0.05)`=plt2,
             `DpnII`=plt.dpnII,
            `T2D Credible Set \n(DIAGRAM)`=plt1,
            `Genes`=ap,
            heights=c(0.80,1.5,0.80,3,1),title=loc) + scale_x_sequnit("Mb") 
    } else{
      tracks(`ATAC \nEndoC`=plt3,
             #`Islet eQTL\n(FDR < 0.05)`=plt2,
             `DpnII`=plt.dpnII,
            `T2D Credible Set \n(DIAGRAM)`=plt1,
            `Genes`=ap,
            heights=c(1,1,3,1),title=loc) + scale_x_sequnit("Mb")    
    }

  } else{
    if(dim(e.sub)[1]>0){
      tracks(`ATAC \nEndoC`=plt3,
             `Islet eQTL\n(FDR < 0.05)`=plt2,
             `DpnII`=plt.dpnII,
             `T2D Credible Set \n(DIAGRAM)`=plt1,
              heights=c(1,1,0.5,2),title=loc) + scale_x_sequnit("Mb") 
    } else{
      tracks(`ATAC \nEndoC`=plt3,
             `DpnII`=plt.dpnII,
             `T2D Credible Set \n(DIAGRAM)`=plt1,
              heights=c(1,0.5,3),title=loc) + scale_x_sequnit("Mb")       
    }
  }
}



