
#' Plot the results from a comparative analysis of 4C/CapC experiment
#'
#' @param x1 data structure containing the 4C/CapC scores (column1: pos, column2: coverage score)
#' @param x2 data structure for second condition (same as x1)
#' @param window number of fragments in a single window
#' @param cut.off quantile cut-off
#' @param abs.cut.off
#' @param vp.pos position of the viewpoint (only necessary if vp.dist is defined)
#' @param vp.dist minimal distance to the viewpoint to consider differential regions
#' @param ... arguments be passed on to the plot command
#'
#' @return data structure with fragments with increased contacts. Column 1: fragment position, column 2: 1 (increased contact in x1) or 2 (increased contact in x2)
#' @export
#'
#' @examples
compare.data.plot <- function(x1, x2, window = 21, cut.off = 0.997, abs.cut.off = 200, vp.pos = 0, vp.dist=0, ...){
  x1 <- x1[x1[,1]%in%x2[,1],]
  x2 <- x2[x2[,1]%in%x1[,1],]
  diff.window <- compare.data(x1, x2, window=window, cut.off=cut.off, abs.cut.off = abs.cut.off, vp.pos = vp.pos, vp.dist = vp.dist)
  x <- rem(x1[,1], window)
  y1 <- running(x1[,2],window)
  y2 <- running(x2[,2],window)
  plot(x, y1, type='l', col="#4A5F70", lwd=3, ..., xlab="chromosomal position", ylab="4C signal", ylim=c(0,3000))
  lines(x,y2, col="#824E4E", lwd=3)
  sel <- x%in%diff.window[,1]
  segments(x[sel],y1[sel],x[sel],y2[sel], col=rgb(1,0,0, 0.2), lwd=2)
  invisible(diff.window)
}

#permute the data between two experiment data stuctures (x1 & x2)
runmean.perm <- function ( x1, x2, k = 21, iter = 1000 ){
  X <- matrix(NA, ncol=iter, nrow=length(x1)-k+1)
  for(i in 1:iter){
    #randomly select a value from one of the two experiments
    X[,i] <- running(ifelse(runif(length(x1)) > 0.5, x1, x2), k)
  }
  return(X)
}

#permute data points for every position between the multiple
#experiments and recalculate the running mean
multiple.mean.perm <- function( data, window = 21, iter = 1000 ){
  store.mat <- matrix( NA, ncol=iter, nrow=nrow(data)-window+1)
  exp1 <- 1:(ncol(data)/2)
  exp2 <- (ncol(data)/2+1):ncol(data)
  for( i in 1:iter){
    cat(i, "\r")
    shuffled.mat <- t(apply(data, 1, function(x) sample(x, length(x)) ))
    store.mat[,i] <- running(apply(shuffled.mat[,exp1],1,median), window)
  }
  store.mat
}

#' Select differentially contacted regions using replicates
#'
#' @param data list containing the coverage scores for the replicates, note that the number of replicate must be the same for both conditions
#' @param window number of fragments in a single window
#' @param cut.off quantile cut-off used for calling a differential site
#' @param abs.cut.off minimal absolute difference between two window to be considered a differential region
#'
#' @return data structure with fragments with increased contacts. Column 1: fragment position, column 2: 1 (increased contact in x1) or 2 (increased contact in x2)
#' @export
#'
#' @examples
compare.data.multiple <- function(data, window = 21, cut.off=0.997, abs.cut.off = 200 ){
  num.exp <- length(data)
  #merge the elements of the data list
  data.m <- data[[1]]
  for( i in 2:num.exp ){
    data.m <- merge(data.m, data[[i]], by=1)
  }
  #store the fragment positions
  pos <- data.m[,1]
  data.m <- data.m[,-1]
  shuffle.profile <- multiple.mean.perm( data.m, window=window, iter=1000)
  exp1 <- 1:(ncol(data.m)/2) #column indexes experiments of the first condition
  exp2 <- (ncol(data.m)/2+1):ncol(data.m) #column indexes experiments of the second condition
  y1 <- running(apply(data.m[,exp1], 1, median), window)
  y2 <- running(apply(data.m[,exp2], 1, median), window)
  up.i   <- which(apply(y1 > shuffle.profile, 1, mean) > cut.off | apply(y2 < shuffle.profile, 1, mean) > cut.off)
  down.i <- which(apply(y2 > shuffle.profile, 1, mean) > cut.off | apply(y1 < shuffle.profile, 1, mean) > cut.off)

  #select the fragments used for
  diff.x <- abs(y1[up.i] - y2[up.i])
  up.i <- up.i[diff.x > abs.cut.off]
  diff.x <- abs(y1[down.i] - y2[down.i])
  down.i <- down.i[diff.x > abs.cut.off]

  up.pos <- rem(pos,window)[up.i]
  down.pos <- rem(pos,window)[down.i]
  if(length(up.pos) > 0){
    up <- data.frame(pos=up.pos, col=1)
  }else{
    up <- data.frame()
  }
  if(length(down.pos) > 0){
    down <- data.frame(pos=down.pos, col=2)
  }else{
    down <- data.frame()
  }
  rbind(up, down)
}



#' Select differentially contacted regions using single experiments
#'
#' @param x1 data structure containing the 4C/CapC scores (column1: pos, column2: coverage score)
#' @param x2 data structure for second condition (same as x1)
#' @param window number of fragments in a single window
#' @param cut.off quantile cut-off
#' @param abs.cut.off
#' @param vp.pos position of the viewpoint (only necessary if vp.dist is defined)
#' @param vp.dist minimal distance to the viewpoint to consider differential regions
#'
#' @return data structure with fragments with increased contacts. Column 1: fragment position, column 2: 1 (increased contact in x1) or 2 (increased contact in x2)
#' @export
#'
#' @examples
compare.data <- function( x1, x2, window = 21, cut.off = 0.997, abs.cut.off = 500, vp.pos = 0, vp.dist=0){
  #check whether x1 and x2 contain the same fragments
  #if not select the overlapping fragments
  x1 <- x1[x1[,1]%in%x2[,1],]
  x2 <- x2[x2[,1]%in%x1[,1],]

  shuffle.profile <- runmean.perm(x1[,2], x2[,2], k = window)
  #this is written somewhat counter-intuitively:
  #if shuffle profile is higher than the running mean of x2 in more then
  #cut.off fraction of the cases, this mean that x1 is up
  #and vice versa
  #up.i <- which(apply(shuffle.profile > running(x2[,2],window), 1, mean) > cut.off)
  #down.i <- which(apply(shuffle.profile < running(x2[,2],window), 1, mean) > cut.off)
  #alternative less counter-intuitive way
  up.i   <- which(apply(running(x1[,2],window) > shuffle.profile, 1, mean) > cut.off | apply(running(x2[,2],window) < shuffle.profile, 1, mean) > cut.off)
  down.i <- which(apply(running(x2[,2],window) > shuffle.profile, 1, mean) > cut.off | apply(running(x1[,2],window) < shuffle.profile, 1, mean) > cut.off)
  x1.run <- running(x1[,2],window)
  x2.run <- running(x2[,2],window)
  diff.x <- abs(x1.run[up.i] - x2.run[up.i])
  up.i <- up.i[diff.x > abs.cut.off]
  diff.x <- abs(x1.run[down.i] - x2.run[down.i])
  down.i <- down.i[diff.x > abs.cut.off]

  up.pos <- rem(x1[,1],window)[up.i]
  down.pos <- rem(x1[,1],window)[down.i]
  if(vp.dist > 0){
    up.pos <- up.pos[up.pos < vp.pos-vp.dist | up.pos > vp.pos + vp.dist]
    down.pos <- down.pos[down.pos < vp.pos-vp.dist | down.pos > vp.pos + vp.dist]
  }
  if(length(up.pos) > 0){
    up <- data.frame(pos=up.pos, col=1)
  }else{
    up <- data.frame()
  }
  if(length(down.pos) > 0){
    down <- data.frame(pos=down.pos, col=2)
  }else{
    down <- data.frame()
  }
  rbind(up, down)
}


