
# Loading necessary R packages
require(rrcov); require(MASS); require(reshape);require(ggplot2)

#####################################################################################################
## function to calculate cov^(-1/2) ##
#####################################################################################################
sqrtMat <- function(x, n){
  with(eigen(x), vectors %*% (values^n * t(vectors)))
}

#####################################################################################################
## Calculating the Robust MCD-based Mahalanobis Distance Outlyingness measure ##
#####################################################################################################
getRMDO <- function(
  dat,
  nsamp="best", # the nsamp variable for CovMcd
  myalpha = 0.5 # alpha for proportion of obs to include in MCD-MD
){
  # Using fast-MCD to calculate robust mean and covariance
  mcd1 <- CovMcd(dat,alpha=myalpha,nsamp=nsamp) 
  x1 <- mcd1@center
  cn <- mcd1@cov
  
  # Calculating the Mahalanobis distance (MD), and then robust MD outlyingness (RMDO) based on robust mean and covariance
  md <- apply(dat,1,function(i){t(i-x1)%*%solve(cn)%*%(i-x1)})
  rmde <- 1 + sqrt(md)
  rmde1 <- 1/rmde
  rmdo <- 1-rmde1 # outlyingness values
  #   
  #   used <- rownames(dat)[mcd1@best] # the subjects who were included in the MCD calculations
  #   
  # Based on 2016 paper: Contributions to Quadratic Form
  smpsd <- sqrt(diag(cn))
  d <- diag(1/smpsd)
  dsd <- d %*% cn %*% d
  
  trans1 <- sweep(dat,2,x1)
  corrmax <- sqrtMat(dsd,-1/2) %*% d
  contrib <- apply(trans1,1,function(i){corrmax %*% as.matrix(i)})
  w2 <- t(contrib**2)
  colnames(w2) <- colnames(dat)
  
  return(list(RMDO=rmdo, mcdCenter=x1, mcdCov=cn, initDat=dat, myalpha=myalpha, W2=w2))
  # mcdUsed=used, DSD=dsd, CorrMax=corrmax, 
}



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