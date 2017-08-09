#functions for the combined experiment analysis

#Two main functions:
#1) the identification of significant peaks
#2) creating domainograms for peaks

#adopted from the script of Eisinga et al. 2013
#based on work by Koziol et al. 2010
righttailgamma = function(r,k,n) 1 - pgamma(-log(r/(n+1)^k),k,scale=1)


#' Read in a list of wig files
#'
#' @param files vector containg paths to wig files
#' @param vp.pos position of the viewpoint
#' @param window flanking sequence to be considered in the analysis
#'
#' @return A list containing two column matrices with the position and the normalized coverage score
#' @export
#'
#' @examples
readMultipleWig <- function( files, vp.pos, window = 700e3 ){
  data.list <- list()
  i <- 1
  for(f in files){
    d <- readqWig(f, vp.pos=vp.pos, window=window)
    data.list[[i]] <- d
    i <- i+1
  }
  data.list
}

#' Read in a list of matrix files
#'
#' @param files vector containg paths to two column matrix files
#' @param vp.pos position of the viewpoint
#' @param window flanking sequence to be considered in the analysis
#'
#' @return A list containing two column matrices with the position and the normalized coverage score
#' @export
#'
#' @examples
readMultiple <- function( files, vp.pos, window = 700e3 ){
  data.list <- list()
  i <- 1
  for(f in files){
    d <- readMatrix(f, vp.pos=vp.pos, window=window)
    data.list[[i]] <- d
    i <- i+1
  }
  data.list
}

#' Read in a multiple experiment file
#'
#' @param file path to the experiment file
#' @param vp.pos position of the viewpoint
#' @param window flanking sequence to be considered in the analysis
#' @param num.exp Number of experiments
#'
#' @return A list containing two column matrices with the position and the normalized coverage score
#' @export
#'
#' @examples
readMultiColumnFile <- function(file, vp.pos, window=700e3, num.exp=4){
  d <- read.delim(file, h=F, stringsAsFactors=F)
  d <- d[,1:(num.exp+1)] #in case there are weird empty columns
  #normalize to 1M reads
  for( i in 1:num.exp){
    num.reads <- sum(d[,i+1], na.rm=T)
    d[,i+1] <- 1e6*d[,i+1]/num.reads
  }
  d <- d[d[,1] > vp.pos-window & d[,1] < vp.pos+window,]
  #make a list out of the matrix to make it compatible with the peak caller
  data.list <- list()
  for( i in 2:(num.exp+1)){
    data.list[[i-1]] <- d[,c(1,i)]
  }
  data.list
}



#' Combined 4C/Capture C analysis
#'
#' @param files a vector containing the path to the 4C/CapC file(s)
#' @param data list containing the 4C/CapC data in two column format
#' @param type file types: wig: "wig" file, "matrix": matrix file (single column), "multi": multi-column data file, "data": data provided by user
#' @param num.exp number of experiments, used in conjuction with "data" and "multi" in type, otherwise it will be inferred from the number of files
#' @param vp.pos viewpoint position, this can be a single value or a two values to analyse a viewpoint region
#' @param genomic.window flanking sequence to be considered in the analysis
#' @param window number of fragments in a window
#' @param alpha false-discovery rate threshold
#' @param alt method of calculating the rank-product p-value if alt is TRUE a combination of the difference and ratio is used, if FALSE only the ratio is used
#' @param plot logical whether to draw a plot of the data (defaults to false)
#' @param y.min,y.max ylim cut-offs for the generated plot
#'
#' @description Wrapper function for identifying interaction peaks above a background distribution. Either a vector of files or a list of dataset can be entered as input. The viewpoint position is given in the vp.pos argument. For single point analysis such as 4C a single point can be given as vp.pos. For Capture-C data a viewpoint region can be given, this region is excluded from subsequent analysis.
#'
#' @return numeric vector containing the significant fragments
#' @export
#'
#' @examples
#' data <- readMultipleWig( files= files, vp.pos=vp.pos, window=genomic.window)
#' db <- combine.experiments(data=data, num.exp=num.exp, vp.pos=vp.pos)
#' combined.p <- combine.pvalue(p=p.val, window=window)
#' sign.fragments <- significant.fragments(p.value = combined.p, pos = db[,1], window = window, sign = alpha)
#'
combined.analyseC <- function( files, data, type = "", num.exp=0, vp.pos, genomic.window = 700e3, window = 21, alpha=0.01, alt=T, plot=F, y.min= 0, y.max=3000){
  if(num.exp==0){
    num.exp <- length(files)
  }
  if(type=="wig"){
    data <- readMultipleWig( files= files, vp.pos=vp.pos, window=genomic.window)
  }else if(type=="matrix"){
    data <- readMultiple( files= files, vp.pos=vp.pos, window=genomic.window)
  }else if(type=="multi"){
    data <- readMultiColumnFile ( file = files, vp.pos = vp.pos, window = genomic.window, num.exp = num.exp)
  }else if(type=="data"){
    data <- data
  }else{
    stop("Provide data in the correct format")
  }
  db <- combine.experiments(data=data, num.exp=num.exp, vp.pos=vp.pos)
  if(alt){
    p.val <- rank.product.p.alt( data = db, num.exp = num.exp )
  }else{
    p.val <- rank.product.p( data = db, num.exp = num.exp )
  }
  #select the significant fragments
  sign.fragments <- significant.fragments(p.value = p.val, pos = db[,1], window = window, sign = alpha)
  #create a plot that highlights the significant fragments
  if(plot){
    plot.combined.data(db=db[,1:(num.exp+1)], window=window, sign.fragments = sign.fragments, y.min = y.min, y.max = y.max)
  }
  sign.fragments
}

#plot the combined 4C/CapC data, peaks
plot.combined.data <- function(db, window, sign.fragments, y.min = 0, y.max = 3000 ){
  x <- rem(db[,1],window)
  y <- running(apply(db[,-1],1,median),window)
  plot(x, y, type='h', col=ifelse(x%in%sign.fragments, 'red', 'grey'), xlab="chromosomal position", ylab="4C signal", ylim=c(y.min,y.max), axes=F )
  axis(2, at=c(0,y.max), las=2)
  at <- seq(200e3*floor(min(db[,1])/200e3), ceiling(max(db[,1])/200e3)*200e3, by=200e3)
  axis(1, at=seq(0,1e9,by=200e3), lab=sprintf("%.1f", seq(0,1e9,by=200e3)/1e6), cex.axis=1.5)
}

#given a list of files draw a domainogram
#' Domainogram analysis
#'
#' @param files a vector containing the path to the 4C/CapC file(s)
#' @param data list containing the 4C/CapC data in two column format
#' @param type file types: wig: "wig" file, "matrix": matrix file (single column), "multi": multi-column data file, "data": data provided by user
#' @param num.exp: number of experiments, used in conjuction with "multi" in type
#' @param vp.pos viewpoint position
#' @param genomic.window flanking sequence to be considered in the analysis
#' @param max.window biggest window to be considered in the analysis
#' @param plot logical, draw a new plot or not
#' @param offset maximum vertical value to be used for the domainogram
#' @param add vertical increase in the position of domainogram
#' @param max.p maximum p-value in the color scale
#' @param min.p minimum p-value in the color scale
#'
#' @return nothing
#' @export
#'
#' @examples
domainogram.analysis <- function( files, data, type="", num.exp=0, vp.pos, genomic.window = 700e3, max.window=201, plot=T, offset=NA, add=0, max.p=10, min.p=0){
  if(num.exp==0){
      num.exp <- length(files)
  }
  if(type=="wig"){
    data <- readMultipleWig( files= files, vp.pos=vp.pos, window=genomic.window)
  }else if(type=="matrix"){
    data <- readMultiple( files= files, vp.pos=vp.pos, window=genomic.window)
  }else if(type=="multi"){
    data <- readMultiColumnFile ( file = files, vp.pos = vp.pos, window = genomic.window, num.exp = num.exp)
  }else if(type=="data"){
    data <- data
  }else{
    stop("Provide data in the correct format")
  }
  db <- combine.experiments(data=data, num.exp=num.exp, vp.pos=vp.pos)
  #calcute the p-values associated with the rank product
  p.val <- rank.product.p.alt( data = db, num.exp = num.exp )
  #use this p-value vector for the visualization in a domainogram
  domainogram.bpspace(raw.p = p.val, position = db[,1], max.window=max.window, plot=plot, offset=offset, add=add, max.p=max.p, min.p=min.p)
}


#' Combine list of experiments into a matrix with the background model
#'
#' @param data list containing the 4C/CapC data in two column format
#' @param num.exp number of experiments
#' @param vp.pos position of the viewpoint
#'
#' @return merged matrix containing: 1. the position of the fragments, 2:(n+1) the data, (n+1):(n+n+1) background models corresponding to the respective datasets
#' @export
#'
#' @examples
combine.experiments <- function( data, num.exp, vp.pos ){
	#create two element vector containing the viewpoint position
	#if only one viewpoint is given
	if(length(vp.pos) == 1){
		vp.pos <- c(vp.pos,vp.pos)
	}
	vp.pos <- sort(vp.pos)
  data.m <- data[[1]]
  for( i in 2:num.exp ){
    data.m <- merge(data.m, data[[i]], by=1)
  }
  #create the background model for the upstream regions
  data.bg <- data.m
  for( i in 1:num.exp ){
    data.bg[data.m[,1] < vp.pos[1],i+1] <- get.background(data.m[data.m[,1] < vp.pos[1],c(1,i+1)], vp.pos[1] )
  }
  #and for the downstream regions
  for( i in 1:num.exp ){
    data.bg[data.m[,1] > vp.pos[2],i+1] <- get.background(data.m[data.m[,1] > vp.pos[2],c(1,i+1)], vp.pos[2] )
  }
	#if two viewpoint fragments are given set the intervening fragments
	#to zero
  #set background to 1 to prevent NaN in the ratio
	if(vp.pos[1] != vp.pos[2]){
		for( i in 1:num.exp){
			data.m[data.m[,1] >= vp.pos[1] & data.m[,1] <= vp.pos[2],i+1] <- 0
			data.bg[data.m[,1] >= vp.pos[1] & data.m[,1] <= vp.pos[2],i] <- 1
		}
	}
  cbind(data.m, data.bg[,-1])
}

#calculate the p-value associated with the rank-product
#the p-value is the probability that there is no enrichment
#over the background model
#' Calculate the rank product p-value
#'
#' @param data data structure generated by combine.experiments
#' @param num.exp the number of experiments in the dataset
#'
#' @return a vector of p-values for every fragment
#' @export
#'
#' @description Ranks are determined based on the ratio of the actual data compared with the background. Datasets are combined by calculating the product of the ranks. The rank product p-value is estimated under the gamma distribution.
#' @seealso rank.product.p
#' @examples
rank.product.p <- function( data, num.exp ){
  ratios <- data[,2:(num.exp+1)]/data[,(2:(num.exp+1))+num.exp]
  rp <- nrow(data)-apply(ratios,2,rank)+1
  rp <- apply(rp,1,prod)
  p <- righttailgamma(rp,num.exp,length(rp))
}

#calculate the rank product in an alternative way that
#includes both the ratio as well as the absolute difference
#between the actual coverage and the background model
#' Calculate the rank product p-value
#'
#' @param data data structure generated by combine.experiments
#' @param num.exp the number of experiments in the dataset
#'
#' @return a vector of p-values for every fragment
#' @export
#'
#' @description This is the alternative version of the rank.product.p function. Ranks are determined based on the difference and the ratio of the data compared to the background.
#' @seealso rank.product.p
#' @examples
rank.product.p.alt <- function( data, num.exp ){
  ratios <- data[,2:(num.exp+1)]/data[,(2:(num.exp+1))+num.exp]
  diff.v <- data[,2:(num.exp+1)]-data[,(2:(num.exp+1))+num.exp]
  rp1 <- nrow(data)-apply(ratios,2,rank)+1
  rp2 <- nrow(data)-apply(diff.v,2,rank)+1
  rp <- apply(rp1 + rp2, 2, rank)
  rp <- apply(rp,1,prod)
  p <- righttailgamma(rp,num.exp,length(rp))
}

#quick way of generating a vector with the required indexes
multi.seq <- function( start, end ){
  x <- rep(start, end-start+1)->x
  df <- diff(x)
  df <- df + 1
  low <- which(df > 1)
  df[low] <- -diff(c(0,low))+1
  add <- c(0,cumsum(df))
  x + add
}


#combine p-values within a given window
#due to RA Fisher (1925)
#' Combine p-values using Fisher's method
#'
#' @param p a vector of p-values
#' @param window number of p-values to combine
#'
#' @return a numerical vector containing the combined p-value for each window (vector is shortened by the window size - 1)
#' @export
#'
#' @examples
combine.pvalue <- function( p, window=21 ){
  chi <- -2*runsum(log(p),window)
  pchisq( chi, window*2, low=F)
}

#return the significant fragments for a given vector
#of p-values
significant.fragments <- function( p.value, pos, window = 21, sign = 0.01 ){
  p.combined <- combine.pvalue( p.value, window=window)
  #correct the nominal p-value for multiple hypothesis testing
  p.combined <- p.adjust(p.combined, method="fdr")
  #determine the significant windows and select the fragments therein
  index.start <- which(p.combined < sign)
  index.end <- index.start + window - 1
  index <- unique(multi.seq(index.start,index.end))
  pos[index]
}


########################################
#BP DOMAINOGRAMS
########################################
#create domainograms based on the raw p-values generated in the rank-product analysis
domainogram.bpspace <- function(raw.p, position, max.window=201, plot=T, offset=NA, add=0, max.p=10, min.p=0){
  if(plot){
    if(is.na(offset)){
      plot(c(min(position),max(position)), c(0,max.window), type='n', axes=F, xlab="", ylab="")
    }else{
      plot(c(0,max(position)), c(0,offset), type='n', axes=F, xlab="", ylab="")
    }
  }
  #loop over a range of windows ending at max.window
  #start with an odd number because we want to have odd numbered windows
  for( i in seq(5,max.window,by=2) ){
    #combine p-values using Fisher's method
    comb.p <- combine.pvalue( raw.p, i);
    #calculate the x1 and x2 positions for the rectangles that we will be drawing
    #x1 is halfway between the fragment and the 5' fragment
    #x2 is halfway between the fragment and the 3' fragment
    plot.pos1 <- (rem(position, i)+head(rem(position,i-2),-2))/2
    plot.pos2 <- (rem(position, i)+tail(rem(position,i-2),-2))/2
    #mat.p[(i-5)/2+1,(i:(i+length(comb.p)-1)-i/2)] <- comb.p }
    p.val <- -log10(comb.p)
    p.val[is.na(p.val)] <- 0 #make sure it doesn't crash when 0/0 happens to be in there
    #shift the p-values to adjust to max.p and min.p
    p.val[p.val > max.p] <- max.p
    p.val <- p.val - min.p
    p.val[p.val < 0] <- 0
    new.max.p <- max.p - min.p
    #calculate a black -> red -> yellow color scheme
    red  <- ifelse(p.val > new.max.p/2, 1, (p.val)/(new.max.p/2))
    green <- ifelse(p.val > new.max.p/2, (p.val-new.max.p/2)/(new.max.p/2), 0)
    col <- rgb(red,green,0)
    rect(plot.pos1, i-1+add, plot.pos2, i+1+add, col=col, border=col)
  }
  at <- seq(200e3*floor(min(position)/200e3), ceiling(max(position)/200e3)*200e3, by=200e3)
  axis(1, at=seq(0,1e9,by=200e3), lab=sprintf("%.1f", seq(0,1e9,by=200e3)/1e6), cex.axis=1.5)

}


