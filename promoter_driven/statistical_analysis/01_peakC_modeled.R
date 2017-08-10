"%&%" <- function(a,b) paste0(a,b)

library("dplyr")
library("data.table")


#serv.dir <- "/Users/jtorres/FUSE2/" # Local 
serv.dir <- "/t1-data/user/hugheslab/jtorres/" # CBRG 
peakc.dir <- serv.dir %&% "analysis/wtsa/software/peakC/" # LOCAL 

source(peakc.dir %&% "util.R")
source(peakc.dir %&% "combined_experiment_analysis.R")
source(peakc.dir %&% "experiment_comparison.R")
source(peakc.dir %&% "single_experiment_analysis.R")


proj.dir <- serv.dir %&% "analysis/wtsa/promoter_driven/statistical_analysis/" 
data.dir <- serv.dir %&% "promoter-driven/"
out.dir <- proj.dir %&% "output_files/"

probe.file <- data.dir %&% "oligo-file.txt"
probe.df <- fread(probe.file)
names(probe.df) <- c("Gene","Chr","Ex.Start","Ex.End","Chrom","Start","End","V8","V9")


rfiles.dir <- data.dir %&% "Stat_Package_CB4/Processed_gfc/R_files/" 
rsuffix <- "_Processed_R_analysis.txt"

args = commandArgs(trailingOnly=TRUE)
#my.win <- 11; my.fdr <- 0.01

if (length(args)!=2) {
  stop("Two arguments must be supplied: window and FDR", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  my.win <- args[1]
  my.fdr <- args[2]
}


# Functions 

find_chrom <- function(genename){
  p.df <- filter(probe.df,Gene==genename)
  chrom <- p.df$Chr
  return(chrom)
}

split_frag <- function(s){
  # example s = "chr10:100074092-100074477"
  vec <- strsplit(s,split=":")[[1]]
  c <- vec[1]
  vec2 <- strsplit(vec[2],split="-")[[1]]
  start <- as.integer(vec2[1])
  stop <- as.integer(vec2[2])
  return(list(c,start,stop))
}

add_frag_cols <- function(df){
  frag.vec <- df$fragment
  chr <- c()
  start <- c()
  end <- c()
  for (i in 1:length(frag.vec)){
    f <- frag.vec[i]
    l <- split_frag(f)
    c <- l[[1]]
    s <- l[[2]]
    e <- l[[3]]
    chr <- append(chr,c)
    start <- append(start,s)
    end <- append(end,e)
  }
  out.df <- cbind(chr,start,end,df)
  return(out.df)
}

read_r_capfile <- function(genename){
  print("\nReading R Capure-C file from Telenius pipeline: " %&% genename)
  chrom <- find_chrom(genename)
  fname <- genename %&% rsuffix #"_Processed_R_analysis.txt"
  rdir <- rfiles.dir
  head <- readLines(rdir%&%fname,n=1)
  head <- strsplit(x=head,split="\t")[[1]]
  head <- c("fragment",head)
  df <- fread(rdir %&% fname)
  names(df) <- head
  print("Subsetting to chromosome: " %&% chrom)
  pat <- chrom %&% ":"
  df <- filter(df,grepl(pat,df$fragment))
  print("Appending fragment information columns...")
  out.df <- add_frag_cols(df)
  print("Sorting data frame")
  out.df <- arrange(out.df,start)
  out.df$chr <- as.character(out.df$chr)
  return(out.df)
}

# t <- suppressWarnings(read_r_capfile("ADCY5")) # test 

make_peakc_mat <- function(cap.df, viewpoint, sampname="conD_1_norm",window=1000e3){
  df <- dplyr::select(cap.df,one_of("start",sampname))
  mat <- as.matrix(df)
  max.count <- max(mat[,2])
  #viewpoint <- #mat[(mat[,2]==max.count),][1]
  #print("Viewpoint: " %&% viewpoint)
  max.pos <- viewpoint + window
  min.pos <- viewpoint - window
  sub.mat <- mat[(mat[,1]<=max.pos & mat[,1]>=min.pos),]
  return(sub.mat)
}


make_expdata_list <- function(cap.df, view.point, window=1000e3, samp.vec){ #=c("condD_1_norm","condD_2_norm","condD_3_norm")){
  capclist <- list()
  i = 0
  for (samp in samp.vec){
    i = i+1
    df <- as.data.frame(make_peakc_mat(cap.df, window=window, viewpoint=view.point,sampname=samp))
    capclist[[i]] <- df
  }
  return(capclist)
}

find_viewpoint <- function(genename){
  p.df <- filter(probe.df,Gene==genename)
  vp <- round(abs(mean(c(p.df$End,p.df$Start))))
  return(vp)
}


# Build result file 

build_result_file <- function(n.exp=3, sample.vec=c("conD_1_norm","conD_2_norm","conD_3_norm")){
  out.df <- c() 
  gene.vec <- probe.df$Gene
  pb <- txtProgressBar(min=0,max=length(gene.vec),style=3)
  for (i in 1:length(gene.vec)){
    setTxtProgressBar(pb,i)
    gene <- gene.vec[i]
    df <- suppressWarnings(read_r_capfile(gene)) 
    chrom <- find_chrom(gene)
    vp <- find_viewpoint(gene)
    data <- make_expdata_list(cap.df=df, view.point=vp, samp.vec=sample.vec, window=1000e3)
    sig.fragments <- combined.analyseC(data=data, type="data", num.exp = n.exp, vp.pos = vp,
                                       window = (as.integer(my.win) + 1), 
                                       plot = FALSE,y.max=400, alpha=my.fdr,alt=FALSE) 
    if (length(sig.fragments)>0){
      stack.df <- filter(df,chr==("chr"%&%chrom),start %in% sig.fragments) %>% 
        dplyr::select(one_of("chr","start","end","conD_mean","DIF_conD_minus_conT_mean"))
      names(stack.df)[c(2,3)] <- c("frag.start","frag.end")
      capture <- rep(gene,dim(stack.df)[1])
      stack.df <- cbind(capture,stack.df)
      stack.df$capture <- as.character(stack.df$capture)
      out.df <- rbind(out.df,stack.df)
    }
  }
  write.table(out.df,file=out.dir%&%"peakC-modeled_fdr"%&%my.fdr%&%"_win"%&%my.win%&%".txt",
              sep="\t",quote=FALSE,row.names=FALSE)
  return(out.df)
}

res.df <- build_result_file()
