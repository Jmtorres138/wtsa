

"%&%" <- function(a,b) paste0(a,b)
library("data.table")
library("tidyverse")
library("peaky")

serv.dir1 <- "/home/jason/science/servers/FUSE5/"

work.dir <- serv.dir1 %&% "projects/wtsa/enhancer_driven/first_round/statistical_anlayses/"
write.dir <- work.dir  %&% "peaky_analyses/analysis_files/"

BTS <- readRDS(file=write.dir%&%"BTS.RDS")
omega_to_use <- readRDS(file=write.dir%&%"omega.to.use.RDS")
idmap.df <- fread(write.dir %&% "capture-bait-id-map.txt")

args = commandArgs(trailingOnly=TRUE)
bait <- args[1]

save.name <- write.dir %&% bait %&% ".pky.txt"

relevant_bait = BTS[baitID==bait]
omega_power = omega_to_use
PKS = peaky(relevant_bait, omega_power, iterations=1e6)
P = interpret_peaky(relevant_bait, PKS, omega_power)

write.table(x=P,file=save.name,sep="\t",quote=F,row.names=F)

#par(mfrow=c(3,1))
#zoom = P[abs(P$dist)<1e6]
#plot(x=zoom$dist, xlab="Distance from bait (bp)",
#     y=zoom$residual, ylab="Adjusted readcount")
#
#plot(x=zoom$dist, xlab="Distance from bait (bp)",
#     y=zoom$beta_mean, ylab="Mean contact strength",
#     col="green")
#
#plot(x=zoom$dist, xlab="Distance from bait (bp)",
#     y=zoom$rjmcmc_pos, ylab="MPPC",
#     col="blue")

