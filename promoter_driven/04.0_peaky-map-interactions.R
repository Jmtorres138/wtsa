"%&%" <- function(a,b) paste0(a,b)
#library("data.table")
#library("tidyverse")
library(peaky)
options(bitmapType='cairo')
work.dir <- "/well/mccarthy/users/jason/projects/wtsa/promoter_driven/"
output.dir <- work.dir %&% "peaky_interactions/"
plot.dir <- output.dir %&% "plots/"
args <- commandArgs(trailingOnly=TRUE)
experiment.name <- args[1] #"promoter-capture"
bait.id <- args[2]
BTS <- readRDS(file=output.dir %&% experiment.name %&% ".BTS.RDS")

relevant_bait = BTS[baitID==bait.id]
zoom = relevant_bait[abs(relevant_bait$dist)<1e6]

png(filename = plot.dir %&% "adj-readcount/chain_2/" %&% experiment.name %&%
      "." %&% bait.id %&% ".adj-readcount.png")
  plot(x=zoom$dist,
     y=zoom$residual,
     xlab="Distance from bait (bp)",
     ylab="Adjusted readcount",
     main=paste("Bait",unique(zoom$baitID)))
dev.off()

#Isolate the highest peak and adjusted readcounts within 40 kbp
peak_center = zoom$dist[which.max(zoom$residual)]
zoom$distance_from_peak = abs(zoom$dist-peak_center)
single_peak = zoom[distance_from_peak<4e5,]

#Define the distance function and find the best fit for omega
dist_exp = function(offset, strength, omega_power){
  omega = 10^omega_power
  strength * exp(-(abs(offset * omega)))
}

## Need to try to fit omega for capture, otherwise will default to -3.8 for omega_power
## see Downes et al. 2019

try_fit_omega <- function(single_peak) {
    omega_power_value <- tryCatch(
        {
            message("Fitting omega for capture...")
            fit_omega = nls(residual~dist_exp(distance_from_peak,strength,
              omega_power),data=single_peak,
              start=list(strength=5, omega_power=-3))
            png(filename = plot.dir %&% "peak-decay/" %&% experiment.name %&%
                  "." %&% bait.id %&% ".peak-decay.png")
              plot(single_peak$distance_from_peak, single_peak$residual,
                    main=paste0("Decay of signal for an isolated peak\nEstimate of omega = 10^",
                               round(coefficients(fit_omega)["omega_power"],3)),
                    xlab="Distance from center of peak (bp)",
                    ylab="Adjusted readcount")
              lines(single_peak$distance_from_peak, fitted.values(fit_omega), col="red", lwd=3)
            dev.off()
            round(coefficients(fit_omega)["omega_power"],3) #omega_power value that will be returned
        },
        error=function(cond) {
            message("Failed to fit omega, will set to -3.8:")
            message("Here's the original error message:")
            message(cond)
            # Return value -3.8 in case of error
            return(-3.8)
        },
        warning=function(cond) {
            message("Omega fit caused a warning")
            message("Here's the original warning message:")
            message(cond)
            # Return value -3.8 in case of warning
            return(-3.8)
        },
        finally={
            message("Completed attempt to fit omega value ")
        }
    )
    return(omega_power_value)
}

omega_power <- -3.8 #try_fit_omega(single_peak) # Manually setting to -3.8 # initially set omega to a fixed value
## omega_power <- try_fit_omega(single_peak)  # defaults to -3.8 if fitting procedure fails


relevant_bait = BTS[baitID==bait.id]

PKS = peaky(relevant_bait, omega_power, iterations=10e6) # crank up to 10 million

P = interpret_peaky(relevant_bait, PKS, omega_power)
P$omega.power <- omega_power
write.table(x=P,file=output.dir %&% experiment.name %&% "/chain_2/" %&%
            experiment.name %&% "." %&% bait.id %&%
            ".peaky-output.txt",sep="\t",quote=F,row.names=F,col.names=T)

#write.table(x=P,file=output.dir %&% experiment.name %&% "/vary_omega/" %&%
#            experiment.name %&% "." %&% bait.id %&%
#            ".peaky-output.txt",sep="\t",quote=F,row.names=F,col.names=T) # this time allowing omega value to vary

png(filename = plot.dir %&% "mppc-plots/chain_2/" %&% experiment.name %&%
      "." %&% bait.id %&% ".mppc.png")
  par(mfrow=c(3,1))
  zoom = P[abs(P$dist)<1e6]
  plot(x=zoom$dist, xlab="Distance from bait (bp)",
      y=zoom$residual, ylab="Adjusted readcount")
  plot(x=zoom$dist, xlab="Distance from bait (bp)",
      y=zoom$beta_mean, ylab="Mean contact strength",
      col="green")
  plot(x=zoom$dist, xlab="Distance from bait (bp)",
      y=zoom$rjmcmc_pos, ylab="MPPC",
      col="blue")
dev.off()
