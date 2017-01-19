#####################################################################################################
## Identifying Cut-off for Outliers with simulation ##
#####################################################################################################
idOutliers <- function(
  obsOV, # output from the getRMDO function
  rmdo_nsamp = 500, # number of iterations - same as nsamp in getRMDO.  "best" would be long but potentially do more than default of 500
  nsamps=250, # number of elements within each sample
  nreps=100, # number of samples to take and average over
  thisquant=0.95 # the quantile to use as cut-offs for outliers
){
  # Pulling necessary output from getRMDO
  mcdMn <- obsOV[["mcdCenter"]]
  mcdCov <- obsOV[["mcdCov"]]
  myalpha <- obsOV[["myalpha"]]
  
  # simulating multivariate normal distributions with the robust mean and covaraince from getRMDO
  outquants <- c()
  for(i in 1:nreps){
    tsamp<- mvrnorm(nsamps,mcdMn,mcdCov)
    trmdo <- getRMDO(tsamp,nsamp=rmdo_nsamp,myalpha=myalpha)   # this is where time is problematic
    outquants[i] <- quantile(trmdo[["RMDO"]],thisquant) # saving the quantile value only
  }
  
  # pulling the mean of all simulated quantiles -- will be the cut-off
  simquant <- mean(outquants)
  
  return(list(cutoff=simquant,cutoff_all=outquants))
}
#####################################################################################################