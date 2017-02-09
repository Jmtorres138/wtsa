"%&%" <- function(a,b) paste0(a,b)

library("data.table")
library("ggplot2")

f <- "/Users/jtorres/FUSE/reference/islet/parker.islet.chromHMM.bed"
df <- read.table(f,header=FALSE,stringsAsFactors=FALSE)

str(df)
