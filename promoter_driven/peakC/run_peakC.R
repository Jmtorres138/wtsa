"%&%" <- function(a,b) paste0(a,b)


library("dplyr")
library("data.table")
library("ggplot2")
library("ggrepel")

#peakc.dir <- "/home/wellgen/jtorres/software/peakC/R/"
peakc.dir <- "/Users/jtorres/Google Drive/Science/Projects/wtsa/promoter_driven/peakC/"
source(peakc.dir %&% "util.R")
source(peakc.dir %&% "combined_experiment_analysis.R")
source(peakc.dir %&% "experiment_comparison.R")
source(peakc.dir %&% "single_experiment_analysis.R")

serv.dir <- "/Users/jtorres/FUSE3/"

proj.dir <- serv.dir %&% "promoter-driven/" #"/t1-data/user/hugheslab/jtorres/promoter-driven/"
plot.dir <- proj.dir %&% "peakC/monotonic/plots/"

probe.file <- proj.dir %&% "oligo-file.txt"
probe.df <- fread(probe.file)
names(probe.df) <- c("Gene","Chr","Ex.Start","Ex.End","Chrom","Start","End","V8","V9")


rfiles.dir <- proj.dir %&% "Stat_Package_CB4/Processed_gfc/R_files/"
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

#print(as.integer(my.win) + 1)

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
	print("Reading R Capure-C file from Telenius pipeline: " %&% genename)
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
	return(out.df)
}

make_peakc_mat <- function(cap.df, viewpoint, sampname="conD_1_norm",window=1000e3){
	df <- dplyr::select(cap.df,one_of("start",sampname))
	mat <- as.matrix(df)
	max.count <- max(mat[,2])
	#viewpoint <- #mat[(mat[,2]==max.count),][1]
	print("Viewpoint: " %&% viewpoint)
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


make_poly_df <- function(df){
	poly.x <- c()
	poly.y <- c()
	for (i in 1:dim(df)[1]){
		x1 <- df$start[i]
		x2 <- df$end[i]
		y <- df$conD_mean[i]
		poly.x <- append(poly.x,c(x1,x1,x2,x2))
		poly.y <- append(poly.y,c(0,y,y,0))
	}
	poly.df <- as.data.frame(cbind(poly.x,poly.y))
	return(poly.df)
}


process_index_snpdf <- function(){
	print("Processing data frame of index SNPs that guided experimental design")
	t1 <- fread(proj.dir %&% "gwas-index.bed")
	t2 <- fread(proj.dir %&% "eQTL-index.bed")
	snp.df <- rbind(t1,t2)
	names(snp.df) <- c("CHR","START","END","ANNOT")
	RSID <- c()
	for (i in 1:dim(snp.df)[1]){
		s <- snp.df$ANNOT[i]
		rs <- strsplit(s,"_")[[1]][1]
		RSID <- append(RSID,rs)
	}
	out.df <- cbind(snp.df,RSID)
	return(out.df)
}


check_sigfrags_snps <- function(df,chrom,sig.fragments){
	# This function is mostly for plotting
	frag.start <- c()
	frag.end <- c()
	snp <- c()
	snp.pos <- c()
	sig.df <- filter(df,start %in% sig.fragments)
	index.df <- process_index_snpdf()
	index.df <- select(index.df,one_of("CHR","START","RSID"))
	index.df <- index.df[!duplicated(index.df$RSID),]
	snp.y <- rep(0,dim(index.df)[1])
	index.df <- cbind(index.df,snp.y)
	index.df <- filter(index.df,CHR==("chr"%&%chrom))
	if (dim(index.df)[1] > 0){
		for (i in 1:dim(index.df)[1]){
			pos <- index.df$START[i]
			sname <- index.df$RSID[i]
			for (e in 1:dim(sig.df)[1]){
				start <- sig.df$start[e]
				end <- sig.df$end[e]
				if (pos >= start & pos <= end){
					frag.start <- append(frag.start,start)
					frag.end <- append(frag.end,end)
					snp <- append(snp,sname)
					snp.pos <- append(snp.pos,pos)
				}
			}
		}
		if (length(snp)>0){
			#snp.y <- rep(0,length(snp))
			out.df <- as.data.frame(cbind(frag.start,frag.end,snp,snp.pos))
			out.df$snp.pos <- as.integer(as.character(out.df$snp.pos))
			#out.df$snp.y <- as.integer(as.character(out.df$snp.y))
			num <- length(unique(out.df$snp.pos))
			print("There are " %&% num %&% " index SNPs that map to significantly interacting fragments")
			out.df$frag.start <- as.integer(as.character(out.df$frag.start))
			out.df$frag.end <- as.integer(as.character(out.df$frag.end))
			return(out.df)
		} else{
			print("No index SNPS map to significant fragments")
			return(NULL)
		}
	} else{
		print("No index SNPs are on the chromosome")
		return(NULL)
	}
}

get_contiguous_sigfrags <- function(overlap.frag, df, sig.fragments){
	# This function is mostly for plotting
	# overlap.frag is two element vector c(start,end)
	vec <- c()
	vec <- append(vec,overlap.frag[1])
	sig.df <- filter(df,start %in% sig.fragments) %>% arrange(start)
	print("Checking upstream...")
	index <- match(overlap.frag[1],df$start)
	status <- TRUE
	count <- 0
	while(status==TRUE){
		index <- index - 1
		up.pos <- df$start[index]
		if (up.pos %in% sig.df$start){
			vec <- append(vec,up.pos)
			count <- count + 1
			print("Upstream match count: " %&% count)
		} else{
			status <- FALSE
		}
	}
	print("Checking downstream...")
	index <- match(overlap.frag[1],df$start)
	status <- TRUE
	count <- 0
	while(status==TRUE){
		index <- index + 1
		down.pos <- df$start[index]
		if (down.pos %in% sig.df$start){
			vec <- append(vec,down.pos)
			count <- count + 1
			print("Downstream match count: " %&% count)
		} else{
			status <- FALSE
		}
	}
	print("Total fragments in contiguous interaction: " %&% length(vec))
	vec <- sort(vec)
	return(vec)
}


get_contiguous_sigfrags_df <- function(overlap, df, sig.fragments){
	# This function is mostly for plotting
	out.df <- c()
	for (i in 1:dim(overlap)[1]){
		overlap.frag <- c(overlap$frag.start[i],overlap$frag.end[i])
		FRAG.START <- get_contiguous_sigfrags(overlap.frag, df, sig.fragments)
		I.NUMBER <- rep(i,length(FRAG.START))
		stack.df <- as.data.frame(cbind(FRAG.START,I.NUMBER))
		out.df <- rbind(out.df,stack.df)
	}
	return(out.df)
}


plot_frag_over_cellmean<- function(genename,chrom,df,sig.fragments,viewpoint,window=100e3){
	xlow <- viewpoint - window
	xhigh <- viewpoint + window
	print("Making polygon data frame for fragments in region")
	full.poly.df <- make_poly_df(df)
	saveRDS(full.poly.df, file=plot.dir %&% "dataframes_for_plots/" %&% genename %&% "_full.poly.df.RDS")
	plt <- ggplot(data=full.poly.df, aes(x=poly.x/1e6,y=poly.y)) +
		geom_polygon(color="gray66",fill="gray74") + theme_bw() +
		coord_cartesian(xlim=c(xlow/1e6,xhigh/1e6)) +
		ggtitle("Locus: " %&% genename) +
		xlab("Position on Chromosome " %&% gsub("chr","",chrom) %&% "(Mb)") +
		ylab("Mapped Reads per Fragment")

	if (length(sig.fragments)>0){
		print("Making polygon data frame for significant fragments in region")
		sig.df <- filter(df,start %in% sig.fragments)
		sig.poly.df <- make_poly_df(sig.df)
		saveRDS(sig.poly.df, file=plot.dir %&% "dataframes_for_plots/" %&% genename %&% "_sig.poly.df.RDS")
		plt <- plt + geom_polygon(data=sig.poly.df,aes(x=poly.x/1e6,y=poly.y),
						color="midnightblue",fill="blue2") +
						geom_hline(yintercept=0,color="gray66",size=0.8)
		overlap <- check_sigfrags_snps(df,chrom,sig.fragments)
		if (!is.null(overlap)){
			contig.df <- get_contiguous_sigfrags_df(overlap,df,sig.fragments)
			if (dim(contig.df)[1] > 0){
				peak.df <- filter(df, start %in% contig.df$FRAG.START)
				peak.poly.df <- make_poly_df(peak.df)
				saveRDS(peak.poly.df, file=plot.dir %&% "dataframes_for_plots/" %&% genename %&% "_peak.poly.df.RDS")

				plt <- plt + geom_polygon(data=peak.poly.df,aes(x=poly.x/1e6,y=poly.y),
							color="firebrick3",fill="firebrick1") +
							geom_hline(yintercept=0,color="gray66",size=0.85)
			}
		}
	}

	index.df <- process_index_snpdf()
	index.df <- select(index.df,one_of("CHR","START","RSID"))
	index.df <- index.df[!duplicated(index.df$RSID),]
	index.df <- filter(index.df,CHR==("chr"%&%chrom),START>=xlow,START<=xhigh)
	snp.y <- rep(0,dim(index.df)[1])
	index.df <- cbind(index.df,snp.y)
	saveRDS(index.df, file=plot.dir %&% "dataframes_for_plots/" %&% genename %&% "_index.df.RDS")

	if (dim(index.df)[1] > 0){
		plt <- plt + geom_point(data=index.df, aes(x=START/1e6,y=snp.y),
			shape=23,size=2,color="black",fill="ghostwhite") +
			geom_label_repel(data=index.df,aes(x=START/1e6,y=snp.y,label=RSID),fill="ghostwhite")
	}

	ggsave(filename=plot.dir%&%"peakC_"%&%genename%&%".pdf",plt,width=8,height=4)
	plt
}


get_contiguous_df <- function(df, sig.fragments){
	print("Determining the number of separate contact regions...")
	# overlap.frag is two element vector c(start,end)
	fvec <- sort(sig.fragments)
	fnum <- c()
	count <- 0
	for (f in fvec){
		index <- match(f,df$start)
		prevf <- df$start[index-1]
		nextf <- df$start[index+1]
		if (prevf %in% sig.fragments){
			fnum <- append(fnum,count)
		} else{
			count <- count + 1
			fnum <- append(fnum,count)
		}
	}
	temp.df <- as.data.frame(cbind(fvec,fnum))
	inter.num <- c()
	inter.start <- c()
	inter.end <- c()
	for (i in 1:length(unique(temp.df$fnum))){
		sub.df <- filter(temp.df,fnum==i)
		sub.f <- sub.df$fvec
		start <- sub.f[1]
		end <- df$end[match(sub.f[dim(sub.df)[1]],df$start)]
		inter.num <- append(inter.num,i)
		inter.start <- append(inter.start,start)
		inter.end <- append(inter.end,end)
	}
	out.df <- as.data.frame(cbind(inter.num,inter.start,inter.end))
	return(out.df)
}


get_contiguous_SNPoverlap_df <- function(contig.df,chrom){
	# contig.df from get_contiguous_df
	print("Determining if separate contact regions contain index snps...")
	inter.num <- c()
	inter.start <- c()
	inter.end <- c()
	snp <- c()
	snp.pos <- c()
	index.df <- process_index_snpdf()
	index.df <- select(index.df,one_of("CHR","START","RSID"))
	index.df <- index.df[!duplicated(index.df$RSID),]
	snp.y <- rep(0,dim(index.df)[1])
	index.df <- cbind(index.df,snp.y)
	index.df <- filter(index.df,CHR==("chr"%&%chrom))
	if (dim(index.df)[1] > 0){
		for (i in 1:dim(index.df)[1]){
			rs <- index.df$SNP[i]
			p <- index.df$START[i]
			for (j in 1:dim(contig.df)[1]){
				s <- contig.df$inter.start[j]
				e <- contig.df$inter.end[j]
				if (p >= s & p <= e){
					inter.num <- c(inter.num, j)
					inter.start <- c(inter.start, s)
					inter.end <- c(inter.end, e)
					snp <- c(snp, rs)
					snp.pos <- c(snp.pos, p)
				}
			}
		}
		out.df <- as.data.frame(cbind(inter.num,inter.start,inter.end,snp,snp.pos))
		out.df$inter.num <- as.integer(as.character(out.df$inter.num))
		out.df$inter.start <- as.integer(as.character(out.df$inter.start))
		out.df$inter.end <- as.integer(as.character(out.df$inter.end))
		out.df$snp.pos <- as.integer(as.character(out.df$snp.pos))
		return(out.df)
	} else{
		return(as.data.frame(c()))
	}
}


#test <- get_contiguous_SNPoverlap_df(get_contiguous_df(df, sig.fragments))


sig_interactions <- function(genename, n.exp=3, sample.vec=c("conD_1_norm","conD_2_norm","conD_3_norm")){
	df <- read_r_capfile(genename)
	chrom <- find_chrom(genename)
	print("Looking up viewpoint for experiment probe")
	vp <- find_viewpoint(genename)
	print("Making input list for peakC...")
	data <- make_expdata_list(cap.df=df, view.point=vp, samp.vec=sample.vec, window=1000e3)
	print("Running peakC to find significant Capture-C interactions...")
	sig.fragments <- combined.analyseC(data=data, type="data", num.exp = n.exp, vp.pos = vp,
		window = (as.integer(my.win) + 1), plot = FALSE,y.max=400, alpha=my.fdr,alt=FALSE) # window = 21 # alt=T
	#overlap <- check_sigfrags_snps(df,chrom,sig.fragments)
	plot_frag_over_cellmean(genename,chrom,df,sig.fragments,viewpoint=vp,window=1000e3) # 500e3
	print("Evaluating significant contacts...")
	contig.df <- get_contiguous_df(df, sig.fragments)
	string.vec <- c()
	if (!is.na(contig.df$inter.start[1])){
		for (i in 1:dim(contig.df)[1]){
			s <- contig.df$inter.start[i]
			e <- contig.df$inter.end[i]
			string <- s %&% "_" %&% e
			string.vec <- append(string.vec,string)
		}
		contacts <- paste(string.vec,collapse=",")
		contact.num <- dim(contig.df)[1]
		print("Looking for SNP overlaps...")
		overlap.df <- get_contiguous_SNPoverlap_df(contig.df,chrom)
		string.vec <- c()
		if (dim(overlap.df)[1] > 0){
			for (i in 1:dim(overlap.df)[1]){
				s <- contig.df$inter.start[i]
				e <- contig.df$inter.end[i]
				string <- s %&% "_" %&% e
				string.vec <- append(string.vec,string)
			}
			overlaps <- paste(string.vec,collapse=",")
			overlap.num <- dim(overlap.df)[1]
		} else{
			overlaps <- NA
			overlap.num <- 0
		}
	} else{
		contacts <- NA
		contact.num <- 0
		overlaps <- NA
		overlap.num <- 0
	}
	#contact.num <- dim(contig.df)[1]
	#overlap.num <- dim(overlap.df)[1]
	return(list(chrom, genename, contact.num, contacts, overlap.num, overlaps))
}

create_summary <- function(){
	#probe.df <- fread(probe.file)
	gene.vec <- probe.df$Gene
	chr <- c()
	gene <- c()
	contacts <- c()
	contact.regions <- c()
	snp.overlaps <- c()
	overlap.regions <- c()
	for (gg in gene.vec){
		g.list <- sig_interactions(gg)
		ch <- g.list[[1]]
		g <- g.list[[2]]
		con <- g.list[[3]]
		conint <- g.list[[4]]
		ov <- g.list[[5]]
		ovint <- g.list[[6]]
		chr <- c(chr,ch)
		gene <- c(gene,g)
		contacts <- c(contacts,con)
		contact.regions <- c(contact.regions,conint)
		snp.overlaps <- c(snp.overlaps,ov)
		overlap.regions <- c(overlap.regions,ovint)
	}
	out.df <- as.data.frame(cbind(chr,gene,contacts,contact.regions,snp.overlaps,overlap.regions))
	out.df$contacts <- as.integer(as.character(out.df$contacts))
	out.df$snp.overlaps <- as.integer(as.character(out.df$snp.overlaps))
	write.table(out.df,file=proj.dir%&%"peakC-summary.window"%&%my.win%&%".fdr"%&%my.fdr%&%".txt",quote=FALSE,sep="\t",row.names=FALSE)
	return(out.df)
}

res.df <- create_summary()
