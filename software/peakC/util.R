
# Helper functions for running the peak calling software.

#running mean function
#' Title
#'
#' @param x numeric vector
#' @param n window size for the running window
#'
#' @return a numeric vector with the windowed means
#' @export
#'
#' @examples
running<-function(x,n=20){
  cumsum(x)->sum.v
  sum.v<-c(0,sum.v)
  #(sum.v[(n+1):length(x)]-sum.v[1:(length(x)-n)])/n
  diff(sum.v,n)/n
}

#running sum
runsum<-function(x,n=20){
  cumsum(x)->sum.v
  sum.v<-c(0,sum.v)
  diff(sum.v,n)
}

#make running function compatible vector
#' Remove leading and trailing values
#'
#' @param a vector
#' @param n window size of the corresponding running function
#' @description remove n/2 elements from the front and the end
#' @return vector, shortened by (n - 1) elements
#' @export
#'
#' @examples
rem <- function(a, n ){
  half.window <- floor(n/2)
  head(tail(a, -half.window),-half.window)
}


#generic function for generating the background model from
#a general dataset
data.background <- function( data, vp.pos, fractile=F ){
  data <- cbind(data,0)
  data[data[,1] < vp.pos,3] <- get.background(data[data[,1] < vp.pos,1:2], vp.pos, fractile=fractile )
  data[data[,1] > vp.pos,3] <- get.background(data[data[,1] > vp.pos,1:2], vp.pos, fractile=fractile )
  data
}

#perform pava regression and return the background regression line
get.background <- function( data, vp.pos, weight.factor=0, fractile=F){
  require(isotone)
  switched = FALSE
  weights <- (1:nrow(data))**weight.factor
  if(data[1,1] > vp.pos){
    data[,1] <- -data[,1] #reverse the sign to make the trend increasing
    switched = TRUE
    weights <- rev(weights)
  }
  #create the isotonic regression
  if(fractile){
    lm <- gpava(data[,1], data[,2], solver=weighted.fractile, weights=NULL, p=0.75)
  }else{
    lm <- gpava(data[,1], data[,2], solver=weighted.mean)
  }
  #try it with the fractile version
  #lm <- gpava(data[,1], data[,2], solver=weighted.fractile, p=0.6)
  if(switched)
    pred.data <- data.frame( -lm$z, lm$x )
  else
    pred.data <- data.frame( lm$z, lm$x )

  pred.data[order(pred.data[,1]),2]
}

#qwigly read a wig file (with only one chromosome)
#' Quickly read and normalize a wig formatted file
#'
#' @param file path to a wiggle file
#' @param window genomic window around the viewpoint to read in
#' @param vp.pos postion of the viewpoint in the genome
#'
#' @return a matrix with two columns, position and score
#' @export
#' @description Wrapper function for reading files that are formatted as wig files. The data is also normalized to 1 million sequencing reads.
#' @examples
#' data <- readqWig(file="alpha.wig", window = 700e3, vp.pos = 32224333 )
readqWig <- function( file, window, vp.pos ){
  wig <- scan(file, skip = 2)
  d <- matrix(wig, ncol=2, byrow=T)
  d <- d[-which.max(d[,2]),]
  d <- d[d[,1]!=vp.pos,]
  #select only the non-blind fragments
  if(sum(d[,2]) > 0){
    d[,2] <- 1e6*d[,2]/sum(d[,2])
  }else{
      stop("Data file does not contain any data")
  }
  if(window > 0){
    d <- d[d[,1] > vp.pos - window & d[,1] < vp.pos + window,]
  }else{
    #select a genomic region around the viewpoint with a given amount of coverage
    i <- range(which( running(d[,2]>0,2001) > 0.2))+1000
    print(i)
    if(any(is.infinite(i))){
      window = 100e3
      d <- d[d[,1] > vp.pos - window & d[,1] < vp.pos + window,]
    }else{
      d <- d[i[1]:i[2],]
    }

  }
  d
}

#qwigly read a 4C file (with only one chromosome)
#' Quickly read and normalize a matrix formatted file
#'
#' @param file path to a two column file
#' @param window genomic window around the viewpoint to read in
#' @param vp.pos postion of the viewpoint in the genome
#'
#' @return a matrix with two columns, position and score
#' @export
#' @description Wrapper function for reading files that are formatted as wig files. The data is also normalized to 1 million sequencing reads.
#' @examples
#'
readMatrix <- function( file, window, vp.pos ){
  vec <- scan(file)
  d <- matrix(vec, ncol=2, byrow=T)
  d <- d[-which.max(d[,2]),]
  d <- d[d[,1]!=vp.pos,]
  #select only the non-blind fragments
  if(sum(d[,2]) > 0){
    d[,2] <- 1e6*d[,2]/sum(d[,2])
  }else{
    stop("Data file does not contain any data")
  }
  if(window > 0){
    d <- d[d[,1] > vp.pos - window & d[,1] < vp.pos + window,]
  }else{
    #select a genomic region around the viewpoint with a given amount of coverage
    i <- range(which( running(d[,2]>0,2001) > 0.2))+1000
    print(i)
    if(any(is.infinite(i))){
      window = 100e3
      d <- d[d[,1] > vp.pos - window & d[,1] < vp.pos + window,]
    }else{
      d <- d[i[1]:i[2],]
    }

  }
  d
}

rem.top.score <- function( x ){
  i <- rev(order(x))
  x[i[1]] <- x[i[2]]
  mean(x)
}

find.consecutive <- function( i, min.diff=10 ){
  start <- c(1,which(diff(i) > min.diff)+1)
  end <- c(which(diff(i) > min.diff), length(i))
  cbind(i[start],i[end])
}

shuffle_along <- function( x, k ){
  x.s <- seq_along(x)
  y <- sample(x.s)
  x[unlist(split(y, (match(y, x.s)-1) %/% k), use.names = FALSE)]
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

