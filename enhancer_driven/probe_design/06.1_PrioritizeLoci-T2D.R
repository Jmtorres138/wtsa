
"%&%" <- function(a,b) paste0(a,b) 
library("data.table")
library("dplyr")

#serv.dir <- "/Users/jtorres/FUSE/"
serv.dir <- "/well/got2d/jason/"

work.dir <- serv.dir %&% "projects/wtsa/enhancer_driven/probe_design/"
rds.dir <- work.dir %&% "rds/"
profile.dir <- work.dir %&% "profile-snps/"

source(serv.dir %&% "projects/wtsa/enhancer_driven/probe_design/" %&%
         "04.2_TrackPlot.R")

save.dir <- work.dir %&% "profile-loci/"
dir.create(save.dir)

fig.dir <- work.dir %&% "figures/"




prof.df <- fread(profile.dir%&%"profile_credt2d.txt")





cumm_99_old <- function(vec){
  # vec is a arranged (descending) vector of PPAs 
  # returns the indeces needed to get cummulative PPA of 0.99 
  sum = 0 
  out.vec <- c()
  for (i in 1:length(vec)){
    val <- vec[i]
    sum <- sum+val
    if (sum<=0.99){
      out.vec<-append(out.vec,i)
    }
  }
  return(out.vec)
}

cumm_99 <- function(vec){
  # vec is a arranged (descending) vector of PPAs 
  # returns the indeces needed to get cummulative PPA of 0.99 
  sum = 0 
  count = 0
  out.vec <- c()
  for (i in 1:length(vec)){
    val <- vec[i]
    sum <- sum+val
    if (sum<=0.99){
      out.vec<-append(out.vec,i)
      count <- count + 1 
    }
  }
  out.vec <- append(out.vec,count+1)
  return(out.vec)
}

build_summary <- function(df=prof.df){
  out.df <- c()
  pb <- txtProgressBar(min=0,max=length(unique(df$SEGNUMBER)),style=3)
  for (i in 1:length(unique(df$SEGNUMBER))){
    setTxtProgressBar(pb,i)
    seg <- unique(df$SEGNUMBER)[i]
    temp <- filter(df,SEGNUMBER==seg) %>% arrange(desc(PPA.fgwas))
    chr <- temp$CHR[1]; locus <- temp$LOCUS[1]; segnum <- temp$SEGNUMBER[1]
    top.endo.atac <- temp$endo.atac[1]; top.num.eqtls <- temp$num.eqtls[1]
    top.egenes <- temp$egenes[1]
    cred.num <- dim(temp)[1]
    fgwas.cred.num <- length(cumm_99(temp$PPA.fgwas))
    if (fgwas.cred.num==0){
      fgwas.cred.num=1
    }
    stack.df <- data.frame(segnum,chr,locus,cred.num,fgwas.cred.num,
                           top.endo.atac,top.num.eqtls,top.egenes)
    out.df <- rbind(out.df,stack.df)
  }
  out.df$chr <- as.character(out.df$chr)
  out.df$locus <- as.character(out.df$locus)
  out.df$top.egenes <- as.character(out.df$top.egenes)
  prop.reduced <- ifelse(out.df$fgwas.cred.num>=1 & out.df$cred.num!=1,
                         (1 - (out.df$fgwas.cred.num/out.df$cred.num)) ,NA)
  out.df <- cbind(out.df,prop.reduced)
  out.df <- arrange(out.df,desc(prop.reduced))
  return(out.df)
}

sum.df <- build_summary()



build_tiers <- function(){
  top.tier <- filter(sum.df,top.endo.atac==TRUE) %>% arrange(desc(top.num.eqtls))
  
  bottom.tier <- filter(sum.df,top.endo.atac==FALSE) %>% 
    arrange(desc(top.num.eqtls))
  
  tier1 <- filter(top.tier,top.num.eqtls>0,fgwas.cred.num<=10)
  tier2 <- filter(top.tier,top.num.eqtls>0,fgwas.cred.num>10)
  tier3 <- filter(top.tier,top.num.eqtls==0,fgwas.cred.num<=10)
  tier4 <- filter(top.tier,top.num.eqtls==0,fgwas.cred.num>10)
  
  tier5 <- filter(bottom.tier,top.num.eqtls>0,fgwas.cred.num<=10)
  tier6 <- filter(bottom.tier,top.num.eqtls>0,fgwas.cred.num>10)
  tier7 <- filter(bottom.tier,top.num.eqtls==0,fgwas.cred.num<=10)
  tier8 <- filter(bottom.tier,top.num.eqtls==0,fgwas.cred.num>10)
  
  check = dim(tier1)[1] +  dim(tier2)[1]+ dim(tier3)[1]+ dim(tier4)[1]+ 
    dim(tier5)[1] + dim(tier6)[1]+ dim(tier7)[1] + 
    dim(tier8)[1] # sanity check   
  if (check==(dim(top.tier)[1]+dim(bottom.tier)[1])){
    tier <- c(rep(1,dim(tier1)[1]), rep(2,dim(tier2)[1]), rep(3,dim(tier3)[1]),
              rep(4,dim(tier4)[1]), rep(5,dim(tier5)[1]), rep(6,dim(tier6)[1]),
              rep(7,dim(tier7)[1]), rep(8,dim(tier8)[1]))
  }
  out.df <- rbind(tier1,tier2,tier3,tier4,tier5,
                  tier6,tier7,tier8)
  out.df <- cbind(out.df,tier)
  write.table(x=out.df,file=save.dir%&%"profile-t2d-loci.txt",
              row.names=F,sep="\t",quote=F)
  return(out.df)
}

tier.df <- build_tiers()

plt <- ggplot(data=filter(tier.df,cred.num<=100),
              aes(x=cred.num,y=fgwas.cred.num)) + 
  geom_point(shape=21,color="black",aes(fill=as.factor(tier))) + theme_bw() + 
  geom_abline(slope=1)


save_top <- function(){
  top.tier <- filter(sum.df,top.endo.atac==TRUE) %>% arrange(desc(top.num.eqtls))
  for (i in 1:length(top.tier$segnum)){
    segnum <- top.tier$segnum[i]
    chr <- top.tier$chr[i]
    locus <- top.tier$locus[i]
    locus <- gsub("/","-",locus); locus <- gsub(" ","-",locus)
    locus <- gsub("(","",locus,fixed = TRUE)
    locus <- gsub(")","",locus,fixed = TRUE)
    savename <- fig.dir %&% segnum %&% "_" %&% chr %&% "_" %&% locus %&% ".png"
    print(paste("index: ",i, " segnum: ",segnum))
    plt <- track_plot(segnum)
    ggsave(savename,plot=p,width=5,height = 8)
  }
}

save_top()




