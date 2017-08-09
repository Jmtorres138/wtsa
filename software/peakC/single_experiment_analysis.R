# functions that perform peak calling on single experiments

#two methods are provided:
# 1) parametric based on the sampled mean
# 2) non-parametric based on the quantile/fractile scores


#take a subset from the data
sample.data <- function( data, fraction = 0.9){
  len <- nrow(data)
  size <- floor(len*fraction)
  i <- sample( 1:len, size)
  data[sort(i),]
}

#find robust scores by random subsampling of the data
sampled.score <- function(d, iter=100, fraction=0.5, window=21, quantile=F){
  sample.mat <- NULL
  #a new matrix is created that contains the sampled running means (#iter) concatenated
  #all in two columns
  for(i in 1:iter){
    d.s <- sample.data(d, fraction=fraction)
    sample.mat <- rbind(sample.mat, cbind(rem(d.s[,1],window), running(d.s[,2],window)))
  }
  #now calculate the sampled score for each position (quantile or min)
  if(quantile)
    sampled.score <- aggregate(sample.mat[,2], list(sample.mat[,1]), quantile, 0.005)
  else
    sampled.score <- aggregate(sample.mat[,2], list(sample.mat[,1]), min)
  sampled.score
}

#spikes are fragments give an outlier result because they show an amplification
#artefact, these are removed from the data set
remove.spikes <- function( s, cut.off=100){
  #find the spikes
  strong.diff <- c(F,head(diff(s[,2]) > cut.off,-1) & tail(diff(s[,2]) < -cut.off, -1),F)
  i <- which(strong.diff)
  #set the value to the average of the flanking fragments
  s[i,2] <- (s[i-1,2]+s[i+1,2])/2
  s
}



#use a non-parametric statistic to identify peaks from 4C data
#the 75% quantile is used as the cut-off for the peaks
fractile.peaks <- function( data, vp.pos, window=21, threshold=2, y.min=500 ){
  data.bg <- data.background( data, vp.pos=vp.pos, fractile=T )
  y.vec <- runquantile(data.bg[,2],k=window,probs=0.75, endrule="constant")

  #add the additional fragments in the window
  index.start <- which(y.vec-data.bg[,3] > y.min & y.vec/data.bg[,3] > threshold) - floor(window/2)
  index.end <- index.start + window - 1

  index <- unique(multi.seq(index.start,index.end))
  index <- index[index > 0]
  #return the fragments that are above the threshold
  data.bg[index,1]
}



#do the pava isotonic regression
#' Perform pava isotonic regression on one side of the data
#'
#' @param data data matrix, containing the up- or downstream 4C/CapC data
#' @param vp.pos positon of the viewpoint in the genome
#' @param window number of fragments to consider in one window
#' @param y.min minimal absolute difference between windowed data and background data
#' @param y.rel minimal relative difference between windowed data and background data
#' @param vp.dist minimal distance a peak needs to be from the viewpoint
#'
#' @return numeric vector containing the fragments in the peaks
#' @export
#'
#' @examples
pava.regression <- function( data, vp.pos, window, y.min = 0, y.rel = 1, vp.dist = 0 ){
  #set the data type to increasing (isotonic) or decreasing (antitonic)
  require(isotone)
  switched = FALSE
  #weights <- (1:nrow(data))**weight.factor
  if(data[1,1] > vp.pos){
    data[,1] <- -data[,1] #reverse the sign to make the trend increasing
    switched = TRUE
    #weights <- rev(weights)
  }
  #create the isotonic regression
  lm <- gpava(data[,1], data[,2], solver=weighted.mean)
  #try it with the fractile version
  #lm <- gpava(data[,1], data[,2], solver=weighted.fractile, p=0.6)
  if(switched)
    pred.data <- data.frame( -lm$z, lm$x )
  else
    pred.data <- data.frame( lm$z, lm$x )

  data[,1] <- abs(data[,1])
  #merge the data to the sampled scores
  sub.score <- merge(data, pred.data, by=1)

  #plot(sub.score[,2], sub.score[,3], pch=20, xlab="actual", ylab="predicted", log='')
  #lm <- lm(sub.score[,3]~
  #lines(lowess(sub.score[,2], sub.score[,3]), col='red')

  #find those regions where the sampled score is above the sampled dataset
  #sig.t <- sub.score[,2] > sub.score[,3] & sub.score[,2]-sub.score[,3] > y.min & sub.score[,3] > 0 & sub.score[,2]/sub.score[,3] > y.rel
  #perform a new selection of sig.t, now y.rel can be < 1 as well, in this way
  #the y.rel becomes a knob that can be turned in such a what that we can get more
  #or less peaks
  #        there should be a difference of y.min    the background model cannot be 0   the ratio between the actual and the background should be over y.rel
  sig.t <- sub.score[,2]-sub.score[,3] > y.min    & sub.score[,3] > 0                  & sub.score[,2]/sub.score[,3] > y.rel

  sig.pos <- sub.score[sig.t,1]

  #get the indexes of the significant fragments
  i <- which(data[,1]%in%sig.pos)
  if(length(i) > 0){
    #collapse them to consecutive fragments
    i.cons <- find.consecutive(i)
    #return the absolute positions
    sig.regions <- cbind(data[i.cons[,1],1], data[i.cons[,2],1])
  }else{
    sig.regions <- NULL
  }

  #select regions that adhere to a minimal set of rules (i.e. distance from the viewpoint and minimal size of 2kb)
  #this prevents that single fragments are scored as interactions
  sig.regions <- sig.regions[abs(sig.regions[,1]-vp.pos) > vp.dist & abs(sig.regions[,2]-vp.pos) > vp.dist & sig.regions[,2]-sig.regions[,1] > 2e3,]
  #fix the problem of the single row matrix matrix becoming a vector
  if(is.null(nrow(sig.regions))){
    if(!is.null(sig.regions)){
      sig.regions <- matrix(sig.regions, ncol=2, byrow=T)
    }else{
      #if sig.regions is NULL turn sig.pos into an empty vector
      sig.regions <- matrix(ncol=2, nrow=0)
    }
  }

  sig.pos <- multi.seq(sig.regions[,1], sig.regions[,2])
  sig.pos <- sig.pos[sig.pos%in%data[,1]]
  #colnames(sig.regions) <- NULL

  sig.pos

}

get.regions <- function( sub.score, y.rel, y.min ){
  #find those regions where the sampled score is above the sampled dataset
  #sig.t <- sub.score[,2] > sub.score[,3] & sub.score[,2]-sub.score[,3] > y.min & sub.score[,3] > 0 & sub.score[,2]/sub.score[,3] > y.rel
  #perform a new selection of sig.t, now y.rel can be < 1 as well, in this way
  #the y.rel becomes a knob that can be turned in such a what that we can get more
  #or less peaks
  #        there should be a difference of y.min    the background model cannot be 0   the ratio between the actual and the background should be over y.rel
  sig.t <- sub.score[,2]-sub.score[,3] > y.min    & sub.score[,3] > 0                  & sub.score[,2]/sub.score[,3] > y.rel
  sig.pos <- sub.score[sig.t,1]
  i <- which(sig.t)

  if(length(i) > 0){
    #collapse them to consecutive fragments
    i.cons <- find.consecutive(i)
    #return the absolute positions
    sig.regions <- cbind(sub.score[i.cons[,1],1], sub.score[i.cons[,2],1])
  }else{
    sig.regions <- NULL
  }

  #select regions that adhere to a minimal set of rules (i.e. distance from the viewpoint and minimal size of 5kb)
  #this prevents that single fragments are scored as interactions
  sig.regions <- sig.regions[sig.regions[,2]-sig.regions[,1] > 5e3,]
  #fix the problem of the single row matrix matrix becoming a vector
  if(is.null(nrow(sig.regions))){
    if(!is.null(sig.regions)){
      sig.regions <- matrix(sig.regions, ncol=2, byrow=T)
    }else{
      #if sig.regions is NULL turn sig.pos into an empty vector
      sig.regions <- matrix(ncol=2, nrow=0)
    }
  }

  sig.pos <- multi.seq(sig.regions[,1], sig.regions[,2])
  #sig.pos <- sig.pos[sig.pos%in%data[,1]]
  colnames(sig.regions) <- NULL

  sig.regions
}

#plot the combined 4C/CapC data, peaks
plot.single.experiment <- function(data, sign.fragments, ... ){
  x <- data[,1]; y <- data[,2]
  plot(x, y, type='h', col=ifelse(x%in%sign.fragments, 'red', 'grey'), xlab="chromosomal position", ylab="4C signal", ylim=c(0,3000), axes=F, ... )
  axis(2, at=c(0,3000), las=2)
  at <- seq(200e3*floor(min(data[,1])/200e3), ceiling(max(data[,1])/200e3)*200e3, by=200e3)
  axis(1, at=at, lab=at/1e6, cex.axis=1.5)
}


#perform the pava regression analysis to identify 4C peaks when starting from a data structure
#' Identify peaks in single 4C experiment
#'
#' @param data a two-column matrix containing the fragment position and the data
#' @param vp.pos viewpoint position
#' @param y.min minimal absolute difference between the background and the experimental data
#' @param y.rel relative difference between the background and the experimental data
#' @param frag.window window size over which the running mean should be calculated
#' @param plot logical whether to draw a plot
#' @param ... addtional arguments passed to the plot that is drawn
#'
#' @return a numeric vector containing the position of the identified peaks
#' @export
#'
#' @examples
single.experiment.peaks <- function( data , vp.pos, y.min = 0, y.rel = 1, frag.window = 21, plot=F, ... ){
  score <- sampled.score(data, iter=1000, fraction=0.8, quantile=T, window=frag.window)
  #remove any spikes that could be in the data
  score <- remove.spikes( score )
  #seperate calculate the data from up- and downstream of the viewpoint
  fragments <- c()
  fragments <- c(fragments, pava.regression(score[score[,1] < vp.pos,], vp.pos, window=frag.window, y.min=y.min, y.rel=y.rel)) #downstream
  fragments <- c(fragments, pava.regression(score[score[,1] > vp.pos,], vp.pos, window=frag.window, y.min=y.min, y.rel=y.rel)) #upstream
  if(plot){
    plot.single.experiment(data=score, sign.fragments = fragments, ... )
  }
  fragments
}
