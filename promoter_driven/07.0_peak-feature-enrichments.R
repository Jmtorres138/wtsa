"%&%" <- function(a,b) paste0(a,b)
library("data.table");library("tidyverse");library("viridis")
library("GenomicRanges")
serv.dir <- "/well/mccarthy/users/jason/"
work.dir <- serv.dir %&% "projects/wtsa/promoter_driven/"
out.dir <- work.dir%&%"output_files/enrichment_tests/"
args = commandArgs(trailingOnly=TRUE)
feature.bed.file <- args[1]
feature.name <- args[2]
feature.df <- fread(feature.bed.file,header=F)
exclude.genes.vec <- c("PAX5","TNFSF11","CR2")
regions.df <- fread(work.dir%&%"output_files/genomic-regions.txt",
                    header=T)
regions.df <- filter(regions.df,!(gene %in% exclude.genes.vec))
## EndoC-betaH1 NG Capture-C peaky interactions (merged)
capc.pky.df <- fread(work.dir%&%"peaky_interactions/"%&%
  "promoter_endo_peaky-interactions.txt")
capc.pky.df$capture <- purrr::map(capc.pky.df$capture,function(s){
  gsub("C15orf38-AP3S2","ARPIN-AP3S2",s)
}) %>% as.character(.)
peaky.df <- capc.pky.df
## Islet pcHi-C peaky interactions (merged)
bed.dir<-work.dir%&%"pcHiC-peaky/"
bed.files <- list.files(bed.dir)
pchic.pky.df <- c()
for (f in bed.files){
  f.df <- fread(bed.dir%&%f)
  s <- strsplit(f,split="_")[[1]][5]
  ##print(s)
  f.df$V4<-s
  names(f.df) <- c("seqnames","start","end","capture")
  pchic.pky.df <- rbind(pchic.pky.df,f.df)
}

message <- function (..., domain = NULL, appendLF = TRUE)
{
    args <- list(...)
    cond <- if (length(args) == 1L && inherits(args[[1L]], "condition")) {
        if (nargs() > 1L)
            warning("additional arguments ignored in message()")
        args[[1L]]
    }
    else {
        msg <- .makeMessage(..., domain = domain, appendLF = appendLF)
        call <- sys.call()
        simpleMessage(msg, call)
    }
    defaultHandler <- function(c) {
        cat(conditionMessage(c), file = stdout(), sep = "")
    }
    withRestarts({
        signalCondition(cond)
        defaultHandler(cond)
    }, muffleMessage = function() NULL)
    invisible()
}

get_gr_obs <- function(capt,feature.df,peaky.df){
  cap.df <- filter(regions.df,gene==capt)
  reg.gr <- GRanges(cap.df$chrom,IRanges(cap.df$start.pos,cap.df$end.pos))
  windows.gr <- slidingWindows(x=reg.gr,width=1000,step=1000)
  feat.df <- filter(feature.df,V1==as.character(seqnames(reg.gr)),
                       V2>=start(reg.gr),V3<=end(reg.gr))
  feature.gr <- GRanges(seqnames=feat.df$V1,
    IRanges(start=feat.df$V2,end=feat.df$V3))
  pky.df <- filter(peaky.df,capture==capt,
    seqnames==as.character(seqnames(reg.gr)),
    start>=start(reg.gr),end<=end(reg.gr))
  pky.gr <- GRanges(seqnames=pky.df$seqnames,
    IRanges(start=pky.df$start,end=pky.df$end))
  return(list(windows.gr,feature.gr,pky.gr))
}
window_cont_table <- function(windows.gr,feature.gr,pky.gr){
  win.gr <- windows.gr[[1]][(width(windows.gr[[1]])>1)]
  cell.a <- 0 # both pky peak and feature
  cell.b <- 0 # No pky peak and feature
  cell.c <- 0 # pky peak but No feature
  cell.d <- 0 # No pky peak and No feature
  ##pb <- txtProgressBar(min=1,max=length(win.gr),style=3)
  for (i in 1:length(win.gr)){
    ##setTxtProgressBar(pb,i)
    in.peak <- win.gr[i] %over% pky.gr
    in.feature <- win.gr[i] %over% feature.gr
    if (in.peak==TRUE&in.feature==TRUE){
      cell.a <- cell.a + 1
    } else if (in.peak==FALSE&in.feature==TRUE){
      cell.b <- cell.b + 1
    } else if (in.peak==TRUE&in.feature==FALSE){
      cell.c <- cell.c + 1
    } else if (in.peak==FALSE&in.feature==FALSE){
      cell.d <- cell.d + 1
    } else{
      stop("Window does not qualify.")
    }
  }
  ##close(pb)
  m <- matrix(c(cell.a,cell.c,cell.b,cell.d),nrow=2,ncol=2)
  rownames(m) <- c("peaky","no.peaky")
  colnames(m) <- c("feature","no.feature")
  #fisher.test(m)
  #chisq.test(m)
  return(m)
}

build_contingency_table <- function(feature.df,feature.name,peaky.df,peaky.name){
  full.mat <- matrix(0,nrow=2,ncol=2)
  for (capt in regions.df$gene){
    if (capt %in% peaky.df$capture){
      message(capt)
      gr.list <- get_gr_obs(capt,feature.df,peaky.df)
      mat <- window_cont_table(gr.list[[1]],gr.list[[2]],gr.list[[3]])
      full.mat <- full.mat + mat
    } else{
      message("peaky results not available for capture: " %&% capt)
    }
  }
  saveRDS(object=full.mat,file=out.dir %&% "contingency-table" %&%
            "_"%&%peaky.name%&%"_"%&%feature.name%&%".RDS")
}
build_contingency_table(feature.df,feature.name,capc.pky.df,"capC")
build_contingency_table(feature.df,feature.name,pchic.pky.df,"pcHiC")
